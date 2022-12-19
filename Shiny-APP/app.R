### 02-05-2022 Update AUC calculation method to default (trapezoid rule)
### 02-05-2022 Update PPFR filter >= 0.5
### 23-05-2022 Update to inlcude the possibility to plot abosulte ratio L/H
### 17-10-2022 V5 Update, include a config.yml file to indicate the ppm mass error, peptide peak found ratio and min number of fragments to use for quant
options(connectionObserver = NULL)
library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(dplyr)
library(ggplot2)
library(ggplot2)
library(data.table)
require(gridExtra)
library(tidyr)
library(ggpubr)
library(MESS)
library(config)

dw <- get("parameter")
ppm = dw$ppm
n_y = dw$n_y
ppfr = dw$ppfr



ui <- fluidPage(

  sidebarPanel(    
    wellPanel(
		tags$h3("Select Dir"),
		shinyDirButton("directory", title=h3("Select Dir"), "Folder select"),verbatimTextOutput("directorypath"),
		tags$hr(),
		textInput("regex", label=p("Does the samples names have a common pattern?")),
		tags$hr(),
		textInput("IT", label=p("Maximum IT used for MSx scans:")),
		tags$hr(),
        fileInput("IS", label = h3("All Precursors mz file")),
        fileInput("t_full", label = h3("Quantification data")),
		uiOutput("samples"),
		uiOutput("sites"),
		tags$hr(),
		#checkboxInput("Relative", "Relative Quantification"),
        checkboxInput("Absolute", "Absolute Quantification"),
		actionButton("goButton", "Plot data selected Phos-site"),
		actionButton("go2Button", "Print Heatmap"),
		downloadButton("downloadData", label = "Download AUC ENDO"),
		downloadButton("downloadData2", label = "Download AUC IS")

	)
  ),
	mainPanel(
			tabsetPanel(type = "tabs",
			tabPanel("Phospho-site Results",
					mainPanel(h3("XIC"), plotOutput("xicPlot", width="100%") %>% withSpinner(color="#0dc5c1")),
					mainPanel(h3("Quantification"), plotOutput("profilePlot") %>% withSpinner(color="#0dc5c1"))
			),
			tabPanel("All results",
				mainPanel(h3("Heatmap"), plotOutput("heatmapPlot", height="900px") %>% withSpinner(color="#0dc5c1"))
			)
			)
	)
)

server <- function(input, output, session) {
		options(shiny.maxRequestSize=90*1024^2)
#FUNCTIONS USED IN THE APP:
#read_sample_table (FUN): reads and process the txt files containing the ITs.		
		read_sample_table<-function(reg_ex,dir_IT, IS, samples){
					p<-read.table(paste(dir_IT,"/",reg_ex,samples,".raw.results.txt",sep=""), head=F, sep="\t")
					p<-p[!duplicated(p),]
					p<-p %>% separate(V3, c("IT.heavy", "IT.light"), ";")
					p$V1<-gsub("\\(", "", p$V1)
					p$V1<-gsub(")", "", p$V1)
					p<-p %>% separate(V1, c("IS", "ENDO"), ",")
					p$IS<-as.numeric(p$IS)
					#Add protein name from IS table loaded in the step before ("All_Precursors_Initial_Survey.csv")
					p<-left_join(p, IS, by="IS")
					colnames(p)<-c("IS", "ENDO", "ScanNumber","IT.heavy", "IT.light", "Site")
					p$IT.heavy<-gsub("IT=", "", p$IT.heavy)
					p$IT.heavy<-as.numeric(p$IT.heavy)
					p$IT.light<-as.numeric(p$IT.light)
					return(p)
		}
#norm_ITs (FUN): normalize the intensity of the heavy and light ions by the injection time. It works for a single fragment at a time (sub_t/sub_p).		
		norm_ITs<-function(sub_t_heavy, sub_t_light, sub_p, IT, filter_l, filter_h) {
					IT<-as.numeric(IT)

					rt_l<-as.double(unlist(strsplit(as.character(sub_t_light$Raw.Times),",")))#rt: vector with retention times
					int_l<-as.double(unlist(strsplit(as.character(sub_t_light$Raw.Intensities),",")))  #int_l: vector with ENDO intensities
					rt_h<-as.double(unlist(strsplit(as.character(sub_t_heavy$Raw.Times),",")))#rt: vector with retention times
					int_h<-as.double(unlist(strsplit(as.character(sub_t_heavy$Raw.Intensities),","))) #int_h: vector with IS intensities
					#filter values that do not match with the ScanIds from p
					rt_l<-rt_l[filter_l]
					int_l<-int_l[filter_l]
					rt_h<-rt_h[filter_h]
					int_h<-int_h[filter_h]
					int_h_norm<-(int_h*(IT/sub_p$IT.heavy)) #int_h_norm: vector with IS intensities corrected by IT
					int_l_norm<-(int_l*(IT/sub_p$IT.light)) #int_l_norm: vector with ENDO intensities corrected by IT		  
					
					ion<-as.character(rep(sub_t_heavy$Fragment.Ion[1], length(rt_h)))
					ion_h<-as.character(rep(sub_t_heavy$Product.Mz[1], length(rt_h)))
					ion_l<-as.character(rep(sub_t_light$Product.Mz[1], length(rt_l)))
					temp<-cbind(rt_l, int_h, int_h_norm, int_l, int_l_norm,ion, ion_h, ion_l)
					return(temp)
		}
#ratio_lh_calculator (FUN): calculates the ratio heavy vs light (for the whole AUC).
		ratio_lh_calculator<-function(sub_t, sub_p, IT) {
			areas_h<-c()
			areas_l<-c()
			IT<-as.numeric(IT)
			scanID<-sub_p$ScanNumber
			sub_t_light<-sub_t[sub_t$Isotope.Label.Type=="light",]
			sub_t_heavy<-sub_t[sub_t$Isotope.Label.Type=="heavy",]
			sub_t_l_ids<-unlist(strsplit(as.character(sub_t_light[1,]$Raw.Spectrum.Ids), ","))
			sub_t_h_ids<-unlist(strsplit(as.character(sub_t_heavy[1,]$Raw.Spectrum.Ids), ","))
			sub_t_l_ids<-gsub(".*\\.", "", sub_t_l_ids)
			sub_t_h_ids<-gsub(".*\\.", "", sub_t_h_ids)
			filter_l<-sub_t_l_ids %in% scanID
			filter_h<-sub_t_h_ids %in% scanID
			for(i in 1:(nrow(sub_t)/2)){
				rt<-c()
				int_h<-c()
				int_l<-c()
				#rt<-as.double(unlist(strsplit(as.character(sub_t$Raw.Times[((nrow(sub_t)/2)+i)]),",")))#rt: vector with retention times
				#int_h<-as.double(unlist(strsplit(as.character(sub_t$Raw.Intensities[((nrow(sub_t)/2)+i)]),","))) #int_h: vector with IS intensities
				#int_l<-as.double(unlist(strsplit(as.character(sub_t$Raw.Intensities[i]),",")))#int_l: vector with ENDO intensities
				rt_l<-as.double(unlist(strsplit(as.character(sub_t_light$Raw.Times[i]),",")))#rt: vector with retention times
				int_l<-as.double(unlist(strsplit(as.character(sub_t_light$Raw.Intensities[i]),",")))  #int_l: vector with ENDO intensities
				rt_h<-as.double(unlist(strsplit(as.character(sub_t_heavy$Raw.Times[i]),",")))#rt: vector with retention times
				int_h<-as.double(unlist(strsplit(as.character(sub_t_heavy$Raw.Intensities[i]),","))) #int_h: vector with IS intensities
				#filter values that do not match with the ScanIds from p
				rt_l<-rt_l[filter_l]
				int_l<-int_l[filter_l]
				rt_h<-rt_h[filter_h]
				int_h<-int_h[filter_h]
				int_h_norm<-(int_h*(IT/sub_p$IT.heavy)) #int_h_norm: vector with IS intensities corrected by IT
				int_l_norm<-(int_l*(IT/sub_p$IT.light)) #int_l_norm: vector with ENDO intensities corrected by IT		  
				areas_h<-c(areas_h, auc(rt_h, int_h_norm)) #areas_h: area of the corrected IS peak for each fragment, and added into a vector that will contain the intensities for all fragment areas.
				areas_l<-c(areas_l, auc(rt_l, int_l_norm)) #areas_l: area of the corrected ENDO peak for each fragment. (...)
			}
			if(length(areas_l)<n_y) {
			  ratio_lh<-NA
			} else {
			  ratio_lh<-sum(areas_l)/sum(areas_h)}
				return(ratio_lh)
		}
#int_calculator (FUN): return the full AUC for the ENDO and IS peptides.		
		int_calculator<-function(sub_t, sub_p, IT) {
			areas_h<-c()
			areas_l<-c()
			IT<-as.numeric(IT)
			scanID<-sub_p$ScanNumber
			sub_t_light<-sub_t[sub_t$Isotope.Label.Type=="light",]
			sub_t_heavy<-sub_t[sub_t$Isotope.Label.Type=="heavy",]
			sub_t_l_ids<-unlist(strsplit(as.character(sub_t_light[1,]$Raw.Spectrum.Ids), ","))
			sub_t_h_ids<-unlist(strsplit(as.character(sub_t_heavy[1,]$Raw.Spectrum.Ids), ","))
			sub_t_l_ids<-gsub(".*\\.", "", sub_t_l_ids)
			sub_t_h_ids<-gsub(".*\\.", "", sub_t_h_ids)
			filter_l<-sub_t_l_ids %in% scanID
			filter_h<-sub_t_h_ids %in% scanID

			for(i in 1:(nrow(sub_t)/2)){
				#rt<-as.double(unlist(strsplit(as.character(sub_t$Raw.Times[((nrow(sub_t)/2)+i)]),",")))#rt: vector with retention times
				#int_h<-as.double(unlist(strsplit(as.character(sub_t$Raw.Intensities[((nrow(sub_t)/2)+i)]),","))) #int_h: vector with IS intensities
				#int_l<-as.double(unlist(strsplit(as.character(sub_t$Raw.Intensities[i]),",")))#int_l: vector with ENDO intensities
				rt_l<-as.double(unlist(strsplit(as.character(sub_t_light$Raw.Times[i]),",")))#rt: vector with retention times
				int_l<-as.double(unlist(strsplit(as.character(sub_t_light$Raw.Intensities[i]),",")))  #int_l: vector with ENDO intensities
				rt_h<-as.double(unlist(strsplit(as.character(sub_t_heavy$Raw.Times[i]),",")))#rt: vector with retention times
				int_h<-as.double(unlist(strsplit(as.character(sub_t_heavy$Raw.Intensities[i]),","))) #int_h: vector with IS intensities
				#filter values that do not match with the ScanIds from p
				rt_l<-rt_l[filter_l]
				int_l<-int_l[filter_l]
				rt_h<-rt_h[filter_h]
				int_h<-int_h[filter_h]
				int_h_norm<-(int_h*(IT/sub_p$IT.heavy)) #int_h_norm: vector with IS intensities corrected by IT
				int_l_norm<-(int_l*(IT/sub_p$IT.light)) #int_l_norm: vector with ENDO intensities corrected by IT		  
				areas_h<-c(areas_h, auc(rt_h, int_h_norm)) #areas_h: area of the corrected IS peak for each fragment, and added into a vector that will contain the intensities for all fragment areas.
				areas_l<-c(areas_l, auc(rt_l, int_l_norm)) #areas_l: area of the corrected ENDO peak for each fragment. (...)
			}
				if(length(areas_l)<n_y){
				  areas_l_all<-NA
				  areas_h_all<-NA
				} else {
			  areas_l_all<-sum(areas_l)
				areas_h_all<-sum(areas_h)
				}
				return(list(areas_l_all, areas_h_all))
		}
#ppm_filter (FUN): remove fragments with mass error > 10 ppm
		ppm_filter<-function(quant) {
			quant$Mass.Error.PPM<-as.numeric(as.character(quant$Mass.Error.PPM))
			quant$Mass.Error.PPM<-abs(quant$Mass.Error.PPM)
			x<-quant[is.na(quant$Mass.Error.PPM),]
			if(nrow(x) !=0){
				for(i in 1:nrow(x)){
				  quant<-quant[!(quant$Protein == x$Protein[i] & quant$Replicate == x$Replicate[i] & quant$Fragment.Ion == x$Fragment.Ion[i]),]
				}
			}
			
			x<-quant[quant$Mass.Error.PPM >ppm,]
			if(nrow(x) !=0){
				for(i in 1:nrow(x)){
				  quant<-quant[!(quant$Protein == x$Protein[i] & quant$Replicate == x$Replicate[i] & quant$Fragment.Ion == x$Fragment.Ion[i]),]
				}
			}
			return(quant)
		}
#XIC_plot (FUN): plot the XIC before and after the IT normalization				
		XIC_plot<- function(reg_ex, dir_IT, samples, site, IS, t_full, IT) {
			# p will contain the information for IT for mz IS and ENDO, Scan number and phospho-site
			p<-read_sample_table(reg_ex, dir_IT, IS, samples)
			# t contains the Raw intensities, Raw elution times and fragment ions (light and heavy) for the corresponding sample
			t<-t_full[t_full$Replicate==samples,]
			t<-t[t$Fragment.Ion!="precursor",]
			t<-t[t$Peptide.Peak.Found.Ratio >=ppfr,]
			t<-ppm_filter(t)
			sub_t<-t[t$Protein==site,] #subset of t with a specific site
			sub_p<-p[p$Site==site,]
			sub_p<-na.omit(sub_p)#subset of p with a specific site
			areas_h<-c()
			areas_l<-c()
			final<-c()
			if(nrow(sub_t)==0 || nrow(sub_p)==0){
				plot_1<-ggplot() + theme_void()+ggtitle(samples)
				plot_2<-ggplot() + theme_void()+ggtitle(" ")
				grid.arrange(plot_1, plot_2, nrow=1)
				stop("No valid fragments for this sample and site")
			}
			scanID<-sub_p$ScanNumber
			sub_t_light<-sub_t[sub_t$Isotope.Label.Type=="light",]
			sub_t_heavy<-sub_t[sub_t$Isotope.Label.Type=="heavy",]
			sub_t_l_ids<-unlist(strsplit(as.character(sub_t_light[1,]$Raw.Spectrum.Ids), ","))
			sub_t_h_ids<-unlist(strsplit(as.character(sub_t_heavy[1,]$Raw.Spectrum.Ids), ","))
			sub_t_l_ids<-gsub(".*\\.", "", sub_t_l_ids)
			sub_t_h_ids<-gsub(".*\\.", "", sub_t_h_ids)
			filter_l<-sub_t_l_ids %in% scanID
			filter_h<-sub_t_h_ids %in% scanID
			final<-c()
			for(i in 1:(nrow(sub_t)/2)){
				sub_t_heavy_1<-sub_t_heavy[i,]
				sub_t_light_1<-sub_t_light[i,]
				temp<-norm_ITs(sub_t_heavy_1, sub_t_light_1, sub_p, IT, filter_l, filter_h)
				final<-rbind(final,temp)
				temp<-c()
			}
			final<-as.data.frame(final)
			final$rt<-as.numeric(as.character(final$rt))
			final$int_h<-as.numeric(as.character(final$int_h))
			final$int_l<-as.numeric(as.character(final$int_l))
			final$int_l_norm<-as.numeric(as.character(final$int_l_norm))
			final$int_h_norm<-as.numeric(as.character(final$int_h_norm))			
			# PLOTTING the XIC before and after normalization
			max_y<-max(final$int_h_norm, final$int_l_norm)+1000
			plot_h_norm<-ggplot(final, aes(x=rt, y=int_h_norm))+
						geom_line(aes(group=ion_h, color=interaction(ion,ion_h,sep="-",lex.order=TRUE)), size=1)+
						scale_y_continuous(labels =scales::scientific)+
						theme_bw()+labs(colour="IS")+
						theme(legend.position="bottom")+
						ggtitle("IS after IT norm")
			plot_l_norm<-ggplot(final, aes(x=rt, y=int_l_norm))+
						geom_line(aes(group=ion_l, color=interaction(ion,ion_l,sep="-",lex.order=TRUE)), size=1)+
						scale_y_continuous(labels =scales::scientific)+
						theme_bw()+labs(colour="ENDO")+
						ggtitle("ENDO after IT norm")+
						theme(legend.position="bottom")
			plot_h<-ggplot(final, aes(x=rt, y=int_h))+
						geom_line(aes(group=ion_h, color=interaction(ion,ion_h,sep="-",lex.order=TRUE)), size=1)+
						scale_y_continuous(labels =scales::scientific)+
						theme_bw()+labs(colour="fragm. ions IS")+
						ggtitle("IS before IT norm")+
						theme(legend.position="bottom")
			plot_l<-ggplot(final, aes(x=rt, y=int_l))+
						geom_line(aes(group=ion_l, color=interaction(ion,ion_l,sep="-",lex.order=TRUE)), size=1)+
						scale_y_continuous(labels =scales::scientific)+
						theme_bw()+labs(colour="fragm. ions ENDO")+
						ggtitle("ENDO before IT norm")+
						theme(legend.position="bottom")

			grid.arrange(plot_h, plot_h_norm, plot_l, plot_l_norm, nrow=2)
		}
#heatmap_plot (FUN): plot the ratio H/L normalizaed accross all provided samples.				
		heatmap_plot<- function(reg_ex, dir_IT, IS, t_full, IT) {
			samples<-levels(as.factor(t_full$Replicate))
			site<-levels(as.factor(t_full$Protein))
			results = data.frame(matrix(vector(), length(site), length(samples)+1,
                dimnames=list(c(), c("Site", samples))),
                stringsAsFactors=F)
			results$Site<-levels(as.factor(t_full$Protein))
			for (k in 1:length(samples)){
				p<-read_sample_table(reg_ex,dir_IT, IS, samples[k])
				# t contains the Raw intensities, Raw elution times and fragment ions (light and heavy) for the corresponding sample
				t<-t_full[t_full$Replicate==samples[k],]
				t<-t[t$Fragment.Ion!="precursor",]
				t<-t[t$Peptide.Peak.Found.Ratio >=ppfr,]
				t<-ppm_filter(t)
				t$Protein<-as.factor(t$Protein)
				sites<-levels(droplevels(t$Protein))
				for (j in 1:length(sites)){
					sub_t<-t[t$Protein==sites[j],] #subset of t with a specific site
					sub_p<-p[p$Site==sites[j],]
					sub_p<-na.omit(sub_p)#subset of p with a specific site
					#remove wrong XIC values from sub_t based on scanIds from sub_p
					if(nrow(sub_t)==0){ratio_lh<-NA}
					else if(nrow(sub_p)==0){ratio_lh<-NA}					
					else{ratio_lh<-ratio_lh_calculator(sub_t, sub_p, IT)}
					results[results$Site==sites[j],k+1]<-ratio_lh
				}
			}
			results[is.na(results)] <- 0
			results_scaled<-results
			results_scaled[,2:(length(samples)+1)]<-t(scale(t(results[,2:(length(samples)+1),])))
			results_scaled<-as.data.table(results_scaled)
			results_scaled[is.na(results_scaled)] <- 0
			results_melt<-melt(results_scaled)
			results_melt$sample<-gsub("_0.05.*", "",results_melt$variable)
			results_melt$replicate<-gsub(".*bef_0", "",results_melt$variable)
			results_melt$gene<-gsub(":.*", "",results_melt$Site)
			results_melt$psite<-gsub(".*:", "",results_melt$Site)
			order<-hclust(dist(apply(results_scaled[,2:ncol(results_scaled)],2, as.numeric)))$order
			results_melt$Site<-factor(x=results_melt$Site, levels=results$Site[order], ordered=T)
			
			heatmap_plot<-ggplot(results_melt, aes(x=variable, y=Site))+
							geom_tile(aes(fill=value))+scale_fill_gradient2(high="red", mid="white", low="blue")+
							theme(text = element_text(size=10),
							axis.text.x = element_text(angle=90, hjust=1))
			heatmap_plot
		}

#heatmap_plot_abs (FUN): plot the ratio H/L accross all provided samples.				
		heatmap_plot_abs<- function(reg_ex, dir_IT, IS, t_full, IT) {
			samples<-levels(as.factor(t_full$Replicate))
			site<-levels(as.factor(t_full$Protein))
			results = data.frame(matrix(vector(), length(site), length(samples)+1,
                dimnames=list(c(), c("Site", samples))),
                stringsAsFactors=F)
			results$Site<-levels(as.factor(t_full$Protein))
			for (k in 1:length(samples)){
				p<-read_sample_table(reg_ex,dir_IT, IS, samples[k])
				# t contains the Raw intensities, Raw elution times and fragment ions (light and heavy) for the corresponding sample
				t<-t_full[t_full$Replicate==samples[k],]
				t<-t[t$Fragment.Ion!="precursor",]
				t<-t[t$Peptide.Peak.Found.Ratio >=ppfr,]
				t<-ppm_filter(t)
				t$Protein<-as.factor(t$Protein)
				sites<-levels(droplevels(t$Protein))
				for (j in 1:length(sites)){
					sub_t<-t[t$Protein==sites[j],] #subset of t with a specific site
					sub_p<-p[p$Site==sites[j],]
					sub_p<-na.omit(sub_p)#subset of p with a specific site
					#remove wrong XIC values from sub_t based on scanIds from sub_p
					if(nrow(sub_t)==0){ratio_lh<-NA}
					else if(nrow(sub_p)==0){ratio_lh<-NA}					
					else{ratio_lh<-ratio_lh_calculator(sub_t, sub_p, IT)}
					results[results$Site==sites[j],k+1]<-ratio_lh
				}
			}
			results[is.na(results)] <- 0
			#results_scaled<-results
			#results_scaled[,2:(length(samples)+1)]<-t(scale(t(results[,2:(length(samples)+1),])))
			#results_scaled<-as.data.table(results_scaled)
			#results_scaled[is.na(results_scaled)] <- 0
			results_melt<-melt(results)
			results_melt$sample<-gsub("_0.05.*", "",results_melt$variable)
			results_melt$replicate<-gsub(".*bef_0", "",results_melt$variable)
			results_melt$gene<-gsub(":.*", "",results_melt$Site)
			results_melt$psite<-gsub(".*:", "",results_melt$Site)
			order<-hclust(dist(apply(results[,2:ncol(results)],2, as.numeric)))$order
			results_melt$Site<-factor(x=results_melt$Site, levels=results$Site[order], ordered=T)
			
			heatmap_plot<-ggplot(results_melt, aes(x=variable, y=Site))+
							geom_tile(aes(fill=value))+scale_fill_gradientn(breaks=c(0,0.025,0.05, 0.075, 0.1), colors=c("white", "coral", "coral4", "darkred", "black"))+
							theme(text = element_text(size=10),
							axis.text.x = element_text(angle=90, hjust=1))
			heatmap_plot
		}
#results_table (FUN): export the AUC as a table.			
		results_table<- function(reg_ex, dir_IT, IS, t_full, IT) {
			samples<-levels(as.factor(t_full$Replicate))
			site<-levels(as.factor(t_full$Protein))
			results = data.frame(matrix(vector(), length(site), length(samples)+1,
                dimnames=list(c(), c("Site", samples))),
                stringsAsFactors=F)
			results$Site<-levels(as.factor(t_full$Protein))
			results_h<-results
			results_l<-results
			for (k in 1:length(samples)){
				p<-read_sample_table(reg_ex,dir_IT, IS, samples[k])
				# t contains the Raw intensities, Raw elution times and fragment ions (light and heavy) for the corresponding sample
				t<-t_full[t_full$Replicate==samples[k],]
				t$Protein<-as.factor(t$Protein)
				t<-t[t$Fragment.Ion!="precursor",]
				t<-t[t$Peptide.Peak.Found.Ratio >=ppfr,]
				t<-ppm_filter(t)
				sites<-levels(droplevels(t$Protein))
				for (j in 1:length(sites)){
					sub_t<-t[t$Protein==sites[j],] #subset of t with a specific site
					sub_p<-p[p$Site==sites[j],]
					sub_p<-na.omit(sub_p)#subset of p with a specific site
					if(nrow(sub_t)==0){auc_l<-NA;auc_h<-NA}
					else if(nrow(sub_p)==0){auc_l<-NA;auc_h<-NA}					
					else{auc_l<-int_calculator(sub_t, sub_p, IT)[[1]]
						auc_h<-int_calculator(sub_t, sub_p, IT)[[2]]}
						results_h[results_h$Site==sites[j],k+1]<-auc_h
						results_l[results_l$Site==sites[j],k+1]<-auc_l
				}
			}
			return(list(results_l, results_h))
		}
#profile_plot (FUN): plot the ratio H/L normalizaed accross all provided samples for a particular site.				
		profile_plot<- function(reg_ex, dir_IT, phos_site, IS, t_full, IT) {
			samples<-levels(as.factor(t_full$Replicate))
			results = data.frame(matrix(vector(), 1, length(samples)+1,
                dimnames=list(c(), c("Site", samples))),
                stringsAsFactors=F)
			results$Site<-phos_site
			for (k in 1:length(samples)){
				p<-read_sample_table(reg_ex,dir_IT, IS, samples[k])
				# t contains the Raw intensities, Raw elution times and fragment ions (light and heavy) for the corresponding sample
				t<-t_full[t_full$Replicate==samples[k],]
				t<-t[t$Fragment.Ion!="precursor",]
				t<-t[t$Peptide.Peak.Found.Ratio >=ppfr,]
				t<-ppm_filter(t)
				sub_t<-t[t$Protein==phos_site,] #subset of t with a specific site
				sub_p<-p[p$Site==phos_site,]
				sub_p<-na.omit(sub_p)#subset of p with a specific site
				if(nrow(sub_t)==0){ratio_lh<-NA}	
				else{ratio_lh<-ratio_lh_calculator(sub_t, sub_p, IT)}
				results[results$Site==phos_site,k+1]<-ratio_lh
			}
			results[is.na(results)] <- 0
			results_scaled<-results
			results_scaled[,2:(length(samples)+1)]<-t(scale(t(results[,2:(length(samples)+1),])))
			results_scaled<-as.data.table(results_scaled)
			results_scaled[is.na(results_scaled)] <- 0
			results_melt<-melt(results_scaled)
			results_melt$sample<-gsub("_0.05.*", "",results_melt$variable)
			results_melt$replicate<-gsub("_EGFR.KRAS_.*", "",results_melt$variable)
			results_melt$replicate<-gsub("_0.*", "",results_melt$variable)
			results_melt$gene<-gsub(":.*", "",results_melt$Site)
			results_melt$psite<-gsub(".*:", "",results_melt$Site)
			#results_melt$replicate<-factor(results_melt$replicate, levels=c("X0h", "X1h", "X3h", "X6h", "X12h", "X24h"))
			line_plot<-ggline(results_melt[results_melt$Site==phos_site,], x="replicate", y="value", add=c("mean_se", "jitter"), color="psite")+
					theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1))
			line_plot
		}
#profile_plot_abs (FUN): plot the ratio H/L.				
		profile_plot_abs<- function(reg_ex, dir_IT, phos_site, IS, t_full, IT) {
			samples<-levels(as.factor(t_full$Replicate))
			results = data.frame(matrix(vector(), 1, length(samples)+1,
                dimnames=list(c(), c("Site", samples))),
                stringsAsFactors=F)
			results$Site<-phos_site
			for (k in 1:length(samples)){
				p<-read_sample_table(reg_ex,dir_IT, IS, samples[k])
				# t contains the Raw intensities, Raw elution times and fragment ions (light and heavy) for the corresponding sample
				t<-t_full[t_full$Replicate==samples[k],]
				t<-t[t$Fragment.Ion!="precursor",]
				t<-t[t$Peptide.Peak.Found.Ratio >=ppfr,]
				t<-ppm_filter(t)
				sub_t<-t[t$Protein==phos_site,] #subset of t with a specific site
				sub_p<-p[p$Site==phos_site,]
				sub_p<-na.omit(sub_p)#subset of p with a specific site
				if(nrow(sub_t)==0){ratio_lh<-NA}	
				else{ratio_lh<-ratio_lh_calculator(sub_t, sub_p, IT)}
				results[results$Site==phos_site,k+1]<-ratio_lh
			}
			results[is.na(results)] <- 0
			#results_scaled<-results
			#results_scaled[,2:(length(samples)+1)]<-t(scale(t(results[,2:(length(samples)+1),])))
			#results_scaled<-as.data.table(results_scaled)
			#results_scaled[is.na(results_scaled)] <- 0
			results_melt<-melt(results)
			results_melt$sample<-gsub("_0.05.*", "",results_melt$variable)
			#results_melt$replicate<-gsub("_EGFR.KRAS_.*", "",results_melt$variable)
			results_melt$replicate<-gsub("_0.*", "",results_melt$variable)
			results_melt$gene<-gsub(":.*", "",results_melt$Site)
			results_melt$psite<-gsub(".*:", "",results_melt$Site)
			#results_melt$replicate<-factor(results_melt$replicate, levels=c("X0h", "X1h", "X3h", "X6h", "X12h", "X24h"))
			line_plot<-ggline(results_melt[results_melt$Site==phos_site,], x="replicate", y="value", add=c("mean_se", "jitter"), color="psite")+
					theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1))
			line_plot
		}

		output$samples<-renderUI({
					inFile <- input$t_full
					if (is.null(inFile))
						return(NULL)
					quant_table<-read.csv(inFile$datapath, head=T, sep=",", colClasses = c(Replicate="character"))
					selectInput("samples", "Select Replicate", choices=levels(as.factor(quant_table$Replicate)))
					})

		output$sites<-renderUI({
					inFile_2 <- input$t_full
					if (is.null(inFile_2))
						return(NULL)
					quant_table<-read.csv(inFile_2$datapath, head=T, sep=",", colClasses = c(Replicate="character"))
					selectInput("sites", "Select Phospho-site", choices=levels(as.factor(quant_table$Protein)))
					})
					

		output$directorypath <- renderText({
			if (is.integer(input$directory)) {
			  cat("No directory has been selected (shinyDirChoose)")
			} else {
			  file_path <- parseDirPath(volumes, input$directory)
			  print(file_path)
			}
		  })

		volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
		shinyDirChoose(input, 
						"directory", 
						roots = volumes, 
						session = session, 
						restrictions = system.file(package = "base"), 
						allowDirCreate = FALSE)

		x <- eventReactive(input$goButton, {
											inFile_3 <- input$t_full
											quant_table<-read.csv(inFile_3$datapath, head=T, sep=",", colClasses = c(Replicate="character"))
											inFile_4 <- input$IS
											precursor_list<-read.csv(inFile_4$datapath, head=T, sep=",")
											site<-input$sites
											sample<-input$samples
											dir_IT<-parseDirPath(volumes, input$directory)
											reg_ex<-input$regex
											IT<-input$IT
											XIC_plot(reg_ex, dir_IT, sample, site, precursor_list,quant_table,IT)
											}
							)
		y <- eventReactive(input$goButton, {
											inFile_3 <- input$t_full
											quant_table<-read.csv(inFile_3$datapath, head=T, sep=",",colClasses = c(Replicate="character"))
											inFile_4 <- input$IS
											precursor_list<-read.csv(inFile_4$datapath, head=T, sep=",")
											site<-input$sites
											sample<-input$samples
											dir_IT<-parseDirPath(volumes, input$directory)
											reg_ex<-input$regex
											IT<-input$IT
											if(input$Absolute){profile_plot_abs(reg_ex, dir_IT, site, precursor_list,quant_table, IT)}
											else {profile_plot(reg_ex, dir_IT, site, precursor_list,quant_table, IT)}
											}
							)
		z <- eventReactive(input$go2Button, {
											inFile_3 <- input$t_full
											quant_table<-read.csv(inFile_3$datapath, head=T, sep=",",colClasses = c(Replicate="character"))
											inFile_4 <- input$IS
											precursor_list<-read.csv(inFile_4$datapath, head=T, sep=",")
											site<-input$sites
											sample<-input$samples
											dir_IT<-parseDirPath(volumes, input$directory)
											reg_ex<-input$regex
											IT<-input$IT
											if(input$Absolute){heatmap_plot_abs(reg_ex, dir_IT, precursor_list,quant_table, IT)}
											else{heatmap_plot(reg_ex, dir_IT, precursor_list,quant_table, IT)}
											}
							)
		data_l <- reactive( {
							inFile_3 <- input$t_full
							quant_table<-read.csv(inFile_3$datapath, head=T, sep=",",colClasses = c(Replicate="character"))
							inFile_4 <- input$IS
							precursor_list<-read.csv(inFile_4$datapath, head=T, sep=",")
							site<-input$sites
							sample<-input$samples
							dir_IT<-parseDirPath(volumes, input$directory)
							reg_ex<-input$regex
							IT<-input$IT
							results_table(reg_ex, dir_IT, precursor_list,quant_table, IT)[[1]]
							}
							)
		data_h <- reactive( {
							inFile_3 <- input$t_full
							quant_table<-read.csv(inFile_3$datapath, head=T, sep=",",colClasses = c(Replicate="character"))
							inFile_4 <- input$IS
							precursor_list<-read.csv(inFile_4$datapath, head=T, sep=",")
							site<-input$sites
							sample<-input$samples
							dir_IT<-parseDirPath(volumes, input$directory)
							reg_ex<-input$regex
							IT<-input$IT
							results_table(reg_ex, dir_IT, precursor_list,quant_table, IT)[[2]]
							}
							)
		output$xicPlot <- renderPlot({x()})
		output$profilePlot <- renderPlot({y()})
		output$heatmapPlot <- renderPlot({z()})
		output$downloadData <- downloadHandler(
				filename = "AUC_ENDO.txt",
				content = function(file){
				  write.table(data_l(), file, row.names=F, sep="\t")
				}
			  )
		output$downloadData2 <- downloadHandler(
				filename = "AUC_IS.txt",
				content = function(file){
				  write.table(data_h(), file, row.names=F, sep="\t")
				}
			  )
}

shinyApp(ui = ui, server = server)
