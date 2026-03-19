
zz=0
for (N in c(20,50)){
  for (n in c(10,20)){
   for (vr_int in c(0.1,0.5)){
      for (mult in sqrt(c(0.01,0.1,0.25,0.5,1))){
          for (rh in c(-0.5,-0.2)){
           zz=zz+1
          cat(
            "set.seed(",zz,")\n",
            "num_cluster =",N,"\n",
            "num_subj=",n,"\n",
            "sd_int=sqrt(",vr_int,")\n",
            "multiplier_int=",mult,"\n",
            "rho=",rh,"\n",
            "source(\"source_sim.R\")"
            ,file=paste0("script_",zz,".R"),sep="")


        }
      }}}}


#to get SEs

zz=0
for (N in c(20,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.1,0.5)){
      for (mult in sqrt(c(0.01,0.1,0.25,0.5,1))){
        for (rh in c(-0.5,-0.2)){
            zz=zz+1
          cat(
            "set.seed(",zz,")\n",
            "num_cluster =",N,"\n",
            "num_subj=",n,"\n",
            "sd_int=sqrt(",vr_int,")\n",
            "multiplier_int=",mult,"\n",
            "rho=",rh,"\n",
            "source(\"source_sim_getSe.R\")"
            ,file=paste0("scriptses_",zz,".R"),sep="")


        }
      }}}}



