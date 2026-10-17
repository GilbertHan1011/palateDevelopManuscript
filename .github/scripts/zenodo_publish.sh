#!/usr/bin/env bash
set -euo pipefail

# Archive a GitHub release snapshot to Zenodo (InvenioRDM REST API).
# Env:
#   ZENODO_TOKEN  Zenodo personal access token (deposit:actions:write)
#   TAG           tag name of the release (e.g. v0.1.0)
#   RECORD_ID     (optional) published record id; empty = create first record
#   ZENODO_API    (optional) API base, default https://sandbox.zenodo.org/api

API="${ZENODO_API:-https://zenodo.org/api}"
: "${ZENODO_TOKEN:?ZENODO_TOKEN is not set}"
: "${TAG:?TAG is not set}"
RECORD_ID="${RECORD_ID:-}"
REPO="${GITHUB_REPOSITORY:?GITHUB_REPOSITORY is not set}"
ZIP_FILE="source-${TAG}.zip"
TODAY=$(date +%F)

curl_json=(-H "Authorization: Bearer ${ZENODO_TOKEN}" -H "Content-Type: application/json")
curl_auth=(-H "Authorization: Bearer ${ZENODO_TOKEN}")

# ---------- 1) Download the official source snapshot of the tag ----------
curl -fsSL -o "${ZIP_FILE}" "https://github.com/${REPO}/archive/refs/tags/${TAG}.zip"
echo ">> Downloaded ${ZIP_FILE} ($(stat -c%s "${ZIP_FILE}") bytes)"

# ---------- 2) Build metadata: zenodo-metadata.json + version/date ----------
BASE_META=$(jq -c . zenodo-metadata.json)
META=$(jq -n --argjson m "${BASE_META}" \
  --arg v "${TAG}" --arg d "${TODAY}" \
  '$m + {version: $v, publication_date: $d}')

# ---------- 3) Create draft: first record, or a new version of RECORD_ID ----------
if [ -z "${RECORD_ID}" ]; then
  echo ">> First deposit: creating a new record"
  BODY=$(jq -n --argjson m "${META}" \
    '{access:{record:"public",files:"public"}, files:{enabled:true}, metadata:$m}')
  DRAFT=$(curl -fsS -X POST "${API}/records" "${curl_json[@]}" -d "${BODY}")
else
  echo ">> Creating a new version draft of record ${RECORD_ID}"
  DRAFT=$(curl -fsS -X POST "${API}/records/${RECORD_ID}/versions" "${curl_json[@]}")
  DRAFT_ID=$(jq -r .id <<<"${DRAFT}")
  # Merge inherited metadata with the fresh version/date
  CUR=$(curl -fsS "${API}/records/${DRAFT_ID}/draft" "${curl_auth[@]}")
  MERGED=$(jq -n --argjson cur "${CUR}" --argjson m "${META}" \
    '{access: $cur.access, metadata: ($cur.metadata * $m)}')
  curl -fsS -X PUT "${API}/records/${DRAFT_ID}/draft" "${curl_json[@]}" -d "${MERGED}" >/dev/null
fi

DRAFT_ID=$(jq -r .id <<<"${DRAFT}")
echo ">> Draft ID: ${DRAFT_ID}"

# ---------- 4) Remove any stale files in the draft, then upload ----------
for f in $(curl -fsS "${API}/records/${DRAFT_ID}/draft/files" "${curl_auth[@]}" \
             | jq -r '.entries[]?.key // empty'); do
  echo ">> Deleting stale file ${f}"
  curl -fsS -X DELETE "${API}/records/${DRAFT_ID}/draft/files/${f}" "${curl_auth[@]}" >/dev/null
done

INIT=$(curl -fsS -X POST "${API}/records/${DRAFT_ID}/draft/files" "${curl_json[@]}" \
  -d "[{\"key\": \"${ZIP_FILE}\"}]")

upload_with_retry() {
  local url="$1" file="$2" attempt
  for attempt in 1 2 3 4 5; do
    if curl -fsS -X PUT "$url" \
      -H "Authorization: Bearer ${ZENODO_TOKEN}" \
      -H "Content-Type: application/octet-stream" \
      --data-binary "@${file}"; then
      return 0
    fi
    echo ">> Upload attempt ${attempt} failed, retrying in 60s..."
    sleep 60
  done
  echo ">> Upload failed after 5 attempts" >&2
  return 1
}

upload_with_retry "$(jq -r '.entries[0].links.content' <<<"${INIT}")" "${ZIP_FILE}"
curl -fsS -X POST "$(jq -r '.entries[0].links.commit' <<<"${INIT}")" "${curl_json[@]}" >/dev/null
echo ">> Uploaded ${ZIP_FILE}"

# ---------- 5) Publish (irreversible!) ----------
PUBLISHED=$(curl -fsS -X POST "${API}/records/${DRAFT_ID}/draft/actions/publish" "${curl_json[@]}")
echo ">> Published:"
jq -r '"   Record ID: \(.id)\n   DOI: \(.pids.doi.identifier // "n/a")\n   URL: \(.links.self_html // "n/a")"' <<<"${PUBLISHED}"
