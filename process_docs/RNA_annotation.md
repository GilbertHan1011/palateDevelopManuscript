# Cluster
Four scRNA samples were integrated with scVI, followed by annotating manually. See the [script](20240429_rna_scVI).
![png](img/anno/output_13_1.png)
![4.30_cellanno](img/anno/4.30_cellanno.png)
Detailed Annotation Evidence
![5.30_annotation.png](img/anno/5.30_annotation.png)
# Marker genes

See this [script](script/20240930_RNA_marker).

![Pasted image 20241203105537.png](img/anno/marker.png)

See this [script](script/20240524_rna_differential_gene_hm).
![5.24_figure3_lh_comm1.png](img/anno/5.24_figure3_lh_comm1.png)

# Marker relationship

K14, K6, K5 are important epithelium markers. I explored the relationship between these marker relationship.
![4.29_coexpression.png](img/anno/4.29_coexpression.png)

This result is derived by this [script](20240429_rna_marker_anno).

# Cluster-sample relationship
![Pasted image 20241203105902.png](img/anno/marker_sample.png)

This script is located here.
# Gene module

I also have run gene module with cNMF. The script is [here](script/20240424_cNMF).
![img](img/anno/output_17_0.png)Then I run GO enrichment with these gene modules.
![NMF_enrich.png](img/anno/NMF_enrich.png)The enrich script is [here](script/20240426_cNMF_enriched).
