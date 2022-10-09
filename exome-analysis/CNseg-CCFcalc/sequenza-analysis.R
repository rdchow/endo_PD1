library(sequenza)
samples = c("KR-CH-106","krch114","NP-19-388","S13-26062","S16-18558","S16-21208","S16-7493","S18-18016","S18-23897","S18-8922","S19-17973","X14-1264","X16-3406","X18-362")

for (c in samples){
    mysample = paste("./seqz-out/",c,".bin50.seqz.gz",sep="")
    myout = paste("./seqz-final/",c,sep="")
    test = sequenza.extract(mysample,verbose=TRUE,parallel=12)
    CP = sequenza.fit(test,female=TRUE,mc.cores=12,cellularity=seq(0.2,1,0.01))
    sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = c, out.dir=myout,female=TRUE)
}