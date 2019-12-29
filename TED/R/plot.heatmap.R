plot.heatmap<-function(dat, pdf.name=NULL , cluster, self=T,show.value=F,metric="euc",
						my_palette= colorRampPalette(c("blue", "white", "red"))(n = 254), title="",symkey=F,symbreaks=F,breaks=NULL){
	
	
	
	if(!is.null(pdf.name)) pdf(paste(pdf.name,".pdf",sep=""),pointsize=8,useDingbats=FALSE )
	
	
	
	if(cluster) {
		if(metric=="euc"){
			hc.cols <- hclust(dist(t(dat)),method="ward.D2")
			hc.rows <- hclust(dist(dat),method="ward.D2")	
		}
		
		if(metric=="is.cor"){
			hc.cols <- hclust(as.dist(1-dat),method="ward.D2")
			hc.rows <- hc.cols
		}
	
		if(metric=="use.cor"){
			hc.cols <- hclust(as.dist(1-cor(dat)),method="ward.D2")
			hc.rows <- hclust(as.dist(1-cor(t(dat))),method="ward.D2")
		}
		
		
		if(is.null(breaks)){
			if(show.value) heatmap.2(dat, 
								col=my_palette, 
								density.info="none", 
								trace="none",
								dendrogram='both', 
								symm=self, 
								symkey= symkey,
								symbreaks=symbreaks, 
								scale="none" , 
								Rowv=as.dendrogram(hc.rows), 
								Colv=as.dendrogram(hc.cols),
								margins=c(10,10),
								cellnote= round(dat,2),
     							notecol="black",
     							na.color=par("bg"),
     							main= title)
     			else heatmap.2(dat, 
						col=my_palette, 
						density.info="none", 
						trace="none",
						dendrogram='both', 
						symm=self, 
						symkey= symkey, 
						symbreaks= symbreaks,
						scale="none" , 
						Rowv=as.dendrogram(hc.rows), 
						Colv=as.dendrogram(hc.cols),
						margins=c(10,10),
						main= title)
		}
		else{
			if(show.value) heatmap.2(dat, 
								col=my_palette, 
								density.info="none", 
								trace="none",
								dendrogram='both', 
								symm=self, 
								symkey= symkey,
								symbreaks=symbreaks, 
								scale="none" , 
								Rowv=as.dendrogram(hc.rows), 
								Colv=as.dendrogram(hc.cols),
								margins=c(10,10),
								cellnote= round(dat,2),
     							notecol="black",
     							na.color=par("bg"),
     							main= title,
     							breaks= breaks)
     			else heatmap.2(dat, 
						col=my_palette, 
						density.info="none", 
						trace="none",
						dendrogram='both', 
						symm=self, 
						symkey= symkey, 
						symbreaks= symbreaks,
						scale="none" , 
						Rowv=as.dendrogram(hc.rows), 
						Colv=as.dendrogram(hc.cols),
						margins=c(10,10),
						main= title,
     						breaks= breaks)
		
		}
	}
	
	else {
		if(is.null(breaks)){
			if(show.value) heatmap.2(dat, 
								col=my_palette, 
								density.info="none", 
								trace="none",
								dendrogram='none', 
								symm=self, 
								symkey= symkey,
								symbreaks= symbreaks, 
								scale="none" , 
								Rowv= FALSE, 
								Colv= FALSE,
								margins=c(10,10),
								cellnote= round(dat,2),
     							notecol="black",
     							na.color=par("bg"),
     							main= title)
     			else heatmap.2(dat, 
						col=my_palette, 
						density.info="none", 
						trace="none",
						dendrogram='none', 
						symm=self, 
						symkey= symkey,
						symbreaks= symbreaks, 
						scale="none" , 
						Rowv= FALSE, 
						Colv= FALSE,
						margins=c(10,10),
						main= title)
		}
		else{
			if(show.value) heatmap.2(dat, 
								col=my_palette, 
								density.info="none", 
								trace="none",
								dendrogram='none', 
								symm=self, 
								symkey= symkey,
								symbreaks= symbreaks, 
								scale="none" , 
								Rowv= FALSE, 
								Colv= FALSE,
								margins=c(10,10),
								cellnote= round(dat,2),
     							notecol="black",
     							na.color=par("bg"),
     							main= title,
     							breaks= breaks)
     			else heatmap.2(dat, 
						col=my_palette, 
						density.info="none", 
						trace="none",
						dendrogram='none', 
						symm=self, 
						symkey= symkey,
						symbreaks= symbreaks, 
						scale="none" , 
						Rowv= FALSE, 
						Colv= FALSE,
						margins=c(10,10),
						main= title,
     						breaks= breaks)
		}
	}
				
	if(!is.null(pdf.name))  dev.off()

}






