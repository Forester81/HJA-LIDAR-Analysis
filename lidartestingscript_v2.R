### HJA Lidar comparison



###### MAIN #####
#environment
memory.limit(size = 35000)
memory.size()
setwd('E:\\HJA_Lidar_Testing\\')
library(raster)
library(data.table)
library(ggplot2)

##function

####shift raster
rastshift = function(brk,xvect,yvect){
  #shifts rasters in x and y using all combinations of shiftvect
  #brk = brk2
  start_time <- Sys.time()
  brk3 = brk
  shiftgrid = expand.grid(xvect,yvect) #list of different coordinate shifts
  #rlist = vector(mode="list",nrow(shiftgrid))
  for (k in c(1:nrow(shiftgrid))){
    r = raster::shift(brk,x=shiftgrid[k,1],y=shiftgrid[k,2])
    r2 = projectRaster(r,brk[[1]],method="bilinear") #has to resample when shifts are in units < 1 m'
    names(r2) = paste(names(brk),shiftgrid[k,1],shiftgrid[k,2],sep="_")
    brk3 = stack(brk3,r2)
    print(shiftgrid[k,1])
    print(shiftgrid[k,2])
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  brk3 = brick(brk3)
  end_time <- Sys.time()
  print(end_time - start_time)
  return(brk3)
}


###raster differencing
# differences rasters and produces pertinent stats
rastdif = function(dat,brk,rltn,var){
  r = dat-brk
  names(r) = names(brk)
  ( x.stats <- t(data.frame(x.mean=cellStats(r, "mean"))) )
  ( x.stats <- rbind(x.stats, t(data.frame(x.sd=cellStats(r, "sd")))) )
  ( x.stats <- rbind(x.stats, t(data.frame(x.quant90=quantile(r, probs = c(0.90), type=7,names = FALSE)))))
  ( x.stats <- rbind(x.stats, t(data.frame(x.quant10=quantile(r, probs = c(0.10), type=7,names = FALSE)))))
  ( x.stats <- rbind(x.stats, t(data.frame(x.quant50=quantile(r, probs = c(0.5), type=7,names = FALSE)))))
  ( x.stats <- rbind(x.stats, t(data.frame(x.MAD=cellStats(abs(r),"mean")))))
  ( x.stats <- rbind(x.stats, t(data.frame(x.RMSD=sqrt(cellStats((r^2),"mean"))))))
  x.stats = data.table(t(x.stats))
  x.stats$variable = names(brk)
  x.stats$difference = var
  x.stats$resolution = rltn
  return(x.stats)
}
# #import data
datmask4 = shapefile('E:\\HJA_Lidar_Testing\\HJA_mask4.shp')
# datmask1 = shapefile('E:\\HJA_Lidar_Testing\\HJA_mask.shp') #entire HJA
# be2008 = raster('E:\\HJA_Lidar_Testing\\hja_2008_be.tif')
# be2014 = raster('E:\\HJA_Lidar_Testing\\hja_2014_be.tif')
# be2016 = raster('E:\\HJA_Lidar_Testing\\hja_2016_be.tif')
# 
# 
# # project data to NAD83 UTM since all data is in different coordinate systems
# # takes several minutes
# template = projectExtent(datmask4,datmask4) #makes the template for projection and cropping
# res(template) = 1
# be2014utm = projectRaster(be2014,template,method="bilinear")
# be2008utm = projectRaster(be2008,template,method="bilinear")
# be2016utm = projectRaster(be2016,template,method="bilinear")
# 
# brk2 = brick(be2008utm,be2014utm,be2016utm) #combine them for easier manipulation
# brk2$hja_2014_be = brk2$hja_2014_be * .3048 #convertmeters
# brk2$hja_2016_be = brk2$hja_2016_be * .3048 #convert meters
# 
# writeRaster(brk2, filename="hja_all_be_clip.tif", format="GTiff", overwrite=TRUE)
brk2 = stack("hja_all_be_clip.tif")
brk2 = brick(brk2)
names(brk2) = c("hja_2008_be","hja_2014_be","hja_2016_be")

#shift rasters
# hypoethesis rasters are shifted horizontally
yvect = c(-1,-.75,-.5,-.25,0)
xvect = c(1,.75,.5,.25,0)
brk3 = rastshift(brk2,xvect,yvect) #shift rasters at all combinations of shiftvec
brk3 = brk3[[4:78]] #drop the first 3 
brknames = names(brk3)
brk16 = brk3[[grep("2016",brknames,value=TRUE)]] 
#new08 = projectRaster(new08,brk2$hja_2008_be) 
brk14 = brk3[[grep("2014",brknames,value=TRUE)]]
brk08 = brk3[[grep("2008",brknames,value=TRUE)]]
remove(brk3)

#difference
#hypothesis: 08 is more different from 14 and 16 than 16 and 14 are from each other

dif08_16 = rastdif(brk2$hja_2008_be,brk16,1,"08-16")
dif08_14 = rastdif(brk2$hja_2008_be,brk14,1,"08-14")
dif16_14 = rastdif(brk2$hja_2016_be,brk14,1,"16-14")


#resample and difference
#resampling to 2 m grid is more effective than shifting 

template2 = projectExtent(datmask4,datmask4) #makes the template for projection and cropping
res(template2) = 2
dif08_16.2 = rastdif(resample(brk2$hja_2008_be,template2),resample(brk16,template2),2,"08-16")
dif08_14.2 = rastdif(resample(brk2$hja_2008_be,template2),resample(brk14,template2),2,"08-14")
dif16_14.2 = rastdif(resample(brk2$hja_2016_be,template2),resample(brk14,template2),2,"16-14")



dt = rbindlist(list(dif08_16,dif08_14,dif16_14))
#dt[,5 := NULL] #drop extra column
test = (strsplit(dt$variable,"_"))
test = array(unlist(test),dim=c(5,75))
dt$shift = paste(test[4,],test[5,],sep=",") #add in column that explicitly identifies the level of shift
dt$shift = gsub("[.]1","-1",dt$shift)
dt$shift = gsub("[.]0","-0",dt$shift)
dt$name = paste(dt$shift) #add in a column for naming -- just in case
dt$name =  factor(dt$name)


##noshift TABLE 1
noshiftdat = subset(dt,shift=="0,0")


##scenario analysis -- resolution and shifting
##looking at SD alone #Figure 2
dt2 = subset(dt,x.sd <0.45) # remove most values to zoom in on data
ggplot(dt2,aes(x=reorder(shift,x.sd,sum),y=x.sd))+
  geom_point()+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  facet_wrap(~difference+resolution,nrow = 3)+
  ylab("Standard Deviation (m)")+
  xlab("XYCoordinate Shift (m)")+
  geom_text(aes(label=ifelse(shift=="0,0",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="0.5,-0.5",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="1,0",as.character(shift),'')),hjust=0.5,vjust=-1)
# looks like 08 should shift - 0.5 m, - 0.5 m because in both the 2014 and 2016 difference this shift and resmapling to 2 m reduced variation the most


##mean
dt3 = subset(dt,x.mean > - 0.2) # remove most values to zoom in on data
ggplot(dt3,aes(x=reorder(shift,x.mean,sum),y=x.mean))+
  geom_point()+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  facet_wrap(~difference+resolution,nrow=3)+
  ylab("Mean Difference (m)")+
  xlab("XYCoordinate Shift (m)")+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=1)+
  geom_text(aes(label=ifelse(shift=="0,0",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="0.5,-0.5",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="1,0",as.character(shift),'')),hjust=0.5,vjust=-1)

##MAD (MAE) #Figure 3
dt4 = subset(dt,x.MAD < 0.3) # remove most values to zoom in on data
ggplot(dt4,aes(x=reorder(shift,x.MAD,sum),y=x.MAD))+
  geom_point()+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  facet_wrap(~difference+resolution,nrow=3)+
  ylab("Mean Absolute Difference (m)")+
  xlab("XYCoordinate Shift (m)")+
  geom_text(aes(label=ifelse(shift=="0,0",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="0.5,-0.5",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="0.75,-0.5",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="1,0",as.character(shift),'')),hjust=0.5,vjust=-1)

#RMSE  #Figure 4
dt5 = subset(dt,x.RMSD < 0.56) # remove most values to zoom in on data
ggplot(dt5,aes(x=reorder(shift,x.RMSD,sum),y=x.RMSD))+
  geom_point()+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  facet_wrap(~difference+resolution,nrow=3)+
  ylab("Root Mean Squared Error (m)")+
  xlab("XYCoordinate Shift (m)")+
  geom_text(aes(label=ifelse(shift=="0,0",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="0.5,-0.5",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="0.75,-0.5",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="1,0",as.character(shift),'')),hjust=0.5,vjust=-1)

#50% - Median #figure 5
dt5 = subset(dt,x.quant50 > -0.20) # remove most values to zoom in on data
ggplot(dt5,aes(x=shift,y=x.quant50))+
  geom_point()+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  facet_wrap(~difference+resolution,nrow=3)+
  ylab("Median (m)")+
  xlab("XYCoordinate Shift (m)")+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=1)+
  geom_text(aes(label=ifelse(shift=="0,0",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="0.5,-0.5",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="0.75,-0.5",as.character(shift),'')),hjust=0.5,vjust=-1)+
  geom_text(aes(label=ifelse(shift=="1,0",as.character(shift),'')),hjust=0.5,vjust=-1)



###does it work
#shift the data; test that shift worked -- NEW EVIDENCE SUGGEST 2014 and 2016 are out of alignment
#shift
new16 = raster::shift(brk2$hja_2016_be,x=0.75,y=-0.5) #inverse of the 'best' solution because 2008 will be changed not 14 or 16
new16 = projectRaster(new16,brk2$hja_2016_be) 
new14 = raster::shift(brk2$hja_2014_be,x=0.75,y=-0.5) #inverse of the 'best' solution because 2008 will be changed not 14 or 16
new14 = projectRaster(new14,brk2$hja_2014_be) 


newdiff16 = rastdif(brk2$hja_2008_be,new16,1,"08_16")
newdiff14 = rastdif(brk2$hja_2008_be,new14,1,"08_14")

#bias calc
dtnew = rbindlist(list(newdiff16,newdiff14))
dtnew$shift = "0.75,-0.5"
dtnew$name = "0.75,-0.5"
dtnew = rbindlist(list(dtnew,subset(noshiftdat,difference %in% c("08-14","08-16") & resolution==1)))

#shifting differences Table 3
subset(dt,difference == "08-14" & name == "0.75,-0.5")[,1:7] - newdiff14[,1:7] #shifted 08
subset(dt,difference == "08-16" & name == "0.75,-0.5")[,1:7] - newdiff16[,1:7] #shifted 08

###BIAS correction 
newdiff16_14 = rastdif(new16,new14,1,"08_16")
newdiff08_14 = rastdif(brk2$hja_2008_be,new14,1,"08_16")

#bias correct --- 2014 as reference and MEAN
bias16 = new16 + newdif16_14$x.mean #2016 bias correction
bias08 = brk2$hja_2008_be + newdif08_14$x.mean #2008 bias correction

#### all changes ###

#how do things stack up
biasdif14_16 = rastdif(brk2$hja_2014_be,bias16,1,"14_16") #mean, MAD, median RMSE are lower
biasdif14_08 = rastdif(brk2$hja_2014_be,bias08,1,"14_08") #mean, median, MAD, and RMSE are lower
biasdif14_08_noshift = rastdif(brk2$hja_2014_be,bias08_noshift,1,"14_08-noshift") #although mean, median, MAD, and rmse are lower, the shifted solution is better!

dtbias = rbindlist(list(biasdif14_16,biasdif14_08))

#bias correcting 2016 and 2008   using MEDIAN -- used in the analysis Tables 5 and 6
#correct bias
bias16.2 = brk2$hja_2016_be + dif14_16$x.quant50
bias08.2 = new08 + dif14_08$x.quant50
bias08_noshift.2 = brk2$hja_2008_be +  dif14_08_old$x.quant50
#difference 2014 from bias corrected 2008 and 2016; also shift corrected 2008
biasdif14_16.2 = rastdif(brk2$hja_2014_be,bias16.2,1,"14_16") #mean, MAD, median RMSE are lower
biasdif14_08.2 = rastdif(brk2$hja_2014_be,bias08.2,1,"14_08") #mean, median, MAD, and RMSE are lower

# #bias correction without 2008 shift correction
# biasdif14_08_noshift.2 = rastdif(brk2$hja_2014_be,bias08_noshift.2,1,"14_08-noshift") #although mean, median, MAD, and rmse are lower, the shifted solution is better!
# biasdif16_08_noshift.2 = rastdif(bias16.2,bias08_noshift.2,1,"16_08-noshift") #although mean, median, MAD, and rmse are lower, the shifted solution is better!
# t1 = dif16_08_old[,1:7] - biasdif16_08_noshift.2[,1:7]
# t2 = dif14_08_old[,1:7] - biasdif14_08_noshift.2[,1:7]


#bias correction and shift correction
t3 = dif14_08_old[,1:7] - biasdif14_08.2[,1:7]
t4 = dif16_08_old[,1:7] - biasdif16_08[,1:7]
t5 = dif14_16_old[,1:7] - biasdif14_16.2[,1:7]

#combine the data into table 5
t.all = rbindlist(list(t1,t2,t3,t4,t5))
t.all$shift = c("no shift","no shift","shift","shift","NA")
t.all$comparison = c("2016 - 2008","2014 - 2008","2014 - 2008","2016 - 2008","2014 - 2016")