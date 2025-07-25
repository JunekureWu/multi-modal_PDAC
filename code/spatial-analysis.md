# Generate the spatial context of sub celltype in PDAC TME
### POSTN+ CAF and SPP1+ Macro
```
markers[["POSTN Fibro"]] <- c("POSTN","COL1A1","DCN")
markers[["SPP1 Macro"]] <- c("SPP1","MARCO","CD68")
set.seed(123)
qvalue <- c(6,6,6,12,9)
ST <- c("PDAC1",'PDAC2','PDAC3','PDAC4','PDAC5')

qvalue <- c(9)
ST<-c("PDAC6" , "PDAC7" , "PDAC8"  ,"PDAC9" , "PDAC10",
      "PDAC11", "PDAC12", "PDAC13" ,"PDAC14", "PDAC15",
      "PDAC16" ,"PDAC17", "PDAC18","PDAC19" ,"PDAC20")
for (i in 1:length(ST)) {
  print(paste("It's caculating ",ST[i]," now!",sep = ""))
  load(paste("analysis/ST_",ST[i],"_Bayes.Rdata",sep = ""))
  sce.enhanced <- spatialEnhance(sce,platform = "Visium",q=qvalue,d=15,nrep = 1000,burn.in = 100,save.chain = F)
  sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                  model="xgboost",
                                  feature_names=purrr::reduce(markers, c),
                                  nrounds=0)
  save(sce,sce.enhanced,file = paste("analysis/ST_",ST[i],"_Bayes_enhanced.Rdata",sep = ""))
  sum_counts <- function(sce, features) {
    if (length(features) > 1) {
      colSums(logcounts(sce)[features, ])
    } else {
      logcounts(sce)[features, ]
    }
  }
  
  spot_expr <- purrr::map(markers, function(xs) sum_counts(sce, xs))
  enhanced_expr <- purrr::map(markers, function(xs) sum_counts(sce.enhanced, xs))
  #
  print(paste("Now, enhance feature of ",ST[i]," is completed!",sep = ""))
  #plot function
  plot_expression <- function(sce, expr, title) {
    featurePlot(sce, expr, color=NA) +ggplot2::scale_fill_gradientn(colours = viridis(50, option = "B"))+
      labs(title=title, fill="Log-normalized\nexpression")
  }
  plot_expression_comparison <- function(cell_type) {
    spot.plot <- plot_expression(sce, 
                                 spot_expr[[cell_type]], 
                                 "Spot")
    enhanced.plot <- plot_expression(sce.enhanced,
                                     enhanced_expr[[cell_type]], 
                                     "Enhanced")
    
    (enhanced.plot) + 
      plot_annotation(title=cell_type,
                      theme=theme(plot.title=element_text(size=18)))
  }

  p4<-plot_expression_comparison("POSTN Fibro")
  p5<-plot_expression_comparison("SPP1 Macro")
  
  tmp <- data.frame(spot_expr$`POSTN Fibro`,spot_expr$`SPP1 Macro`)
  colnames(tmp) <-c('POSTN Fibro','SPP1 Macro')
  ggscatter( tmp ,x = 'POSTN Fibro', y = 'SPP1 Macro',
                  add = "reg.line", conf.int = TRUE,color = col1[1],
                  add.params = list(fill = col1[20]),
                  ggtheme = theme_minimal()
  )+ stat_cor(method = "spearman",
              color='black',
              p.accuracy = 0.001
  )

}

```
