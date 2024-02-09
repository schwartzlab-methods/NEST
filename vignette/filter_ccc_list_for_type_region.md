
This workflow will demonstrate how to filter predicted CCC based on cell type, component/region, and ligand-receptor pair. We will start with the visualization without any particular filter and then we will show how each filter works, step by step. 


````
nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name 'NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=3000
````
![png file of the generated altair plot for top 40000 CCC](../images/altair_plot_human_lymph_top3000.png)
![screenshot of the generated histogram plot for top 40000 CCC](../images/histogram_human_lymph_top3000.png)


## Filter by component/region
````
nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name 'NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=3000 --filter=1 --filter_by_component=7
````
![png file of the generated altair plot for top 40000 CCC](../images/altair_plot_human_lymph_top3000_comp7.png)
![screenshot of the generated histogram plot for top 40000 CCC](../images/histogram_human_lymph_top3000_comp7.png)


## Filter by annotation
````
nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=400000 --filter=1 --filter_by_annotation='T-cell' --annotation_file_path=data/V1_Human_Lymph_Node_spatial_annotation.csv 
````
![png file of the generated altair plot for top 40000 CCC](../images/altair_plot_human_lymph_top400000_tcell.png)
![screenshot of the generated histogram plot for top 40000 CCC](../images/histogram_human_lymph_top400000_tcell.png)

## Filter by ligand-receptor pair

```
nest visualize --data_name='V1_Human_Lymph_Node_spatial' --model_name='NEST_V1_Human_Lymph_Node_spatial' --top_edge_count=400000 --filter=1 --filter_by_ligand_receptor='CCL19-CCR7'
```

![png file of the generated altair plot for top 40000 CCC](../images/altair_plot_human_lymph_top400000_ccl19_ccr7.png)
![screenshot of the generated histogram plot for top 40000 CCC](../images/histogram_human_lymph_top400000_ccl19_ccr7.png)
