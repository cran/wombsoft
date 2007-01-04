Ranks<-function(data)
  {
  coord<-as.matrix(data[[1]])
  grid<-data[[3]]
  nb_x<-length(grid[[1]])
  nb_y<-length(grid[[2]])
  ranks<-c()
  for (j in 1:nb_y)
    {
    for (i in 1:nb_x)
      { 
        x_loc<-coord[,1]-grid[[1]][i]
        y_loc<-coord[,2]-grid[[2]][j]
        dist<-x_loc*x_loc+y_loc*y_loc
        rank(dist)->rk
        ranks<-c(ranks,rk)
      }
    }
  ranks
  }

WomblingDominantRk<-function(data,ranks,k,output_syst="syst.txt",output_dir="dir.txt")
  {
  coord<-as.matrix(data[[1]])
  genetic<-as.matrix(data[[2]])
  grid<-data[[3]]
  cvx_mat<-data[[5]]
  nb_ind<-data[[6]]
  nb_indiv_mir<-nrow(coord)
  nb_data<-ncol(genetic)
  nb_x<-length(grid[[1]])
  nb_y<-length(grid[[2]])
  syst<-matrix(0,nb_x,nb_y)
  dire<-matrix(0,nb_x,nb_y)
 .C("wombling_aflp_rg",as.double(ranks),as.double(genetic),as.double(coord),as.integer(nb_indiv_mir),as.integer(nb_data),as.double(grid[[1]]),as.integer(nb_x),as.double(grid[[2]]),as.integer(nb_y),as.double(cvx_mat),as.double(k),as.double(syst),as.double(dire),PACKAGE="wombsoft")->res
  systt<-res[[12]]
  dirr<-res[[13]]
  systemic<-matrix(systt,nb_x,nb_y)
  direction<-matrix(dirr,nb_x,nb_y)
  write.table(systemic,output_syst)
  write.table(direction,output_dir)
  image(grid[[1]],grid[[2]],systemic,zlim=c(min(systt[systt!=0]),max(systt)),col=terrain.colors(200),xlab='x_coordinates',ylab='y_coordinates')
  title(main="Systemic Function")
  points(coord[1:nb_ind,])
  area<-sqrt((grid[[1]][nb_x]-grid[[1]][1])^2+(grid[[2]][nb_y]-grid[[2]][1])^2)
  list_seg<-list(c(),c(),c(),c())
	for (i in 1:((nb_x)/5))
		{
		for (j in 1:((nb_y)/5))
			{
			z<-PlotDir(direction[5*i,5*j],grid[[1]][5*i],grid[[2]][5*j],(4*area/(nb_y))*(systemic[5*i,5*j])/max(systt))
			for (k in 1:4)
				{
				list_seg[[k]]<-c(list_seg[[k]],z[k])
				}
			}
		}
	segments(list_seg[[1]],list_seg[[2]],list_seg[[3]],list_seg[[4]],col='black')
  }
  
 
 
WomblingCodominantRk<-function(data,ranks,k,output_syst="syst.txt",output_dir="dir.txt")
  {
  coord<-as.matrix(data[[1]])
  genetic<-as.matrix(data[[2]])
  grid<-data[[3]]
  cvx_mat<-data[[5]]
  nb_ind_mir<-nrow(coord)         # nombre d'individus distincts !
  nb_indiv<-nrow(genetic)
  nb_data<-2*ncol(genetic)
  nb_x<-length(grid[[1]])
  nb_y<-length(grid[[2]])
  syst<-matrix(0,nb_x,nb_y)
  dire<-matrix(0,nb_x,nb_y)
  nb_all_col<-rep(0,nb_data/2);
  for (i in 1:(nb_data/2))
    {
    nb_all_col[i]<-max(genetic[,i])-min(genetic[,i][genetic[,i]>0])+1
    }
  .C("wombling_codo_rg",as.double(ranks),as.double(genetic),as.double(coord),as.integer(nb_ind_mir),as.integer(nb_data),as.integer(nb_all_col),as.double(grid[[1]]),as.integer(nb_x), as.double(grid[[2]]), as.integer(nb_y),as.double(cvx_mat),as.double(k),as.double(syst),as.double(dire),PACKAGE="wombsoft")->res
  systt<-res[[13]]
  dirr<-res[[14]]  
  systemic<-matrix(systt,nb_x,nb_y)
  direction<-matrix(dirr,nb_x,nb_y)
  write.table(systemic,output_syst)
  write.table(direction,output_dir)
  
  image(grid[[1]],grid[[2]],systemic,zlim=c(min(systt[systt!=0]),max(systt)),col=terrain.colors(200),xlab='x_coordinates',ylab='y_coordinates')
  title(main="Systemic Function")
  points(coord[1:data[[6]],])
  area<-sqrt((grid[[1]][nb_x]-grid[[1]][1])^2+(grid[[2]][nb_y]-grid[[2]][1])^2)
  list_seg<-list(c(),c(),c(),c())
	for (i in 1:((nb_x)/5))
		{
		for (j in 1:((nb_y)/5))
			{
			z<-PlotDir(direction[5*i,5*j],grid[[1]][5*i],grid[[2]][5*j],(4*area/(nb_y))*(systemic[5*i,5*j])/max(systt))
			for (k in 1:4)
				{
				list_seg[[k]]<-c(list_seg[[k]],z[k])
				}
			}
		}
	segments(list_seg[[1]],list_seg[[2]],list_seg[[3]],list_seg[[4]],col='black')
 } 



 CandidateBoundariesCodominantRk<-function(data,ranks,k,pB=0.3)
  {
  coord<-as.matrix(data[[1]])
  genetic<-data[[2]]
  grid<-data[[3]]
  grid_x<-grid[[1]]
  grid_y<-grid[[2]]
  cvx_mat<-data[[5]]
  nb_ind<-data[[6]]              # individus de départ
  nb_ind_mir<-nrow(coord)        # nombre d'individus distincts en comptant ceux recopiés dans le miroir !
  nb_indiv<-nrow(genetic)
  nb_data<-2*ncol(genetic)
  nb_x<-length(grid_x)
  nb_y<-length(grid_y)
  ntot<-nb_x*nb_y
  syst<-matrix(0,nb_x,nb_y)
  dire<-matrix(0,nb_x,nb_y)
  bounds<-matrix(0,nb_x,nb_y)
  dir_bounds<-matrix(0,nb_x,nb_y)
  nb_all_col<-rep(0,nb_data/2)
  nb_alleles<-max(genetic)
  ou_alleles<-c()
  for (i in 1:(nb_data/2))
    {
    mini<-min(genetic[,i][genetic[,i]>0])
    maxi<-max(genetic[,i])
    nb_all_col[i]<-maxi-mini+1
    ou_alleles<-c(ou_alleles,rep(i,nb_all_col[i]))      
    }
  for (k in 1:nb_alleles)
      {
      m<-matrix(0,nb_x,nb_y)
      n<-matrix(0,nb_x,nb_y)
      .C("systemic_allelic_codom_rg",as.double(m),as.double(n),as.double(coord),as.double(ranks),as.double(genetic),as.double(k),as.integer(k),as.integer(ou_alleles[k]),as.double(grid_x),as.double(grid_y),as.integer(nb_x),as.integer(nb_y),as.integer(nb_ind_mir),as.double(cvx_mat),PACKAGE="wombsoft")->res
      syst<-matrix(res[[1]],nb_x,nb_y)
      dire<-matrix(res[[2]],nb_x,nb_y)
      n<-length(syst[syst==0])
		  l_quant<-quantile(syst,probs=c(n/ntot+(1-pB)*(1-n/ntot)))
      for (i in 1:nb_x)
		    {
		    for (j in 1:nb_y)
			   {
			   if (syst[i,j]>l_quant)
			     {
			     bounds[i,j]<-bounds[i,j]+1
			     dir_bounds[i,j]<-dir_bounds[i,j]+2*dire[i,j]
			     }
	       }
        }
        rm(m)
        rm(n)
        }
  list(bounds,dir_bounds,pB)
  }
  
  
 CandidateBoundariesDominantRk<-function(data,ranks,k,pB=0.3)
  {
  coord<-as.matrix(data[[1]])
  genetic<-as.matrix(data[[2]])
  grid<-data[[3]]
  grid_x<-grid[[1]]
  grid_y<-grid[[2]]
  cvx_mat<-data[[5]]
  nb_ind<-data[[6]]              # individus de départ
  nb_ind_mir<-nrow(coord)        # nombre d'individus distincts en comptant ceux recopiés dans le miroir !
  nb_data<-ncol(genetic)
  nb_x<-length(grid_x)
  nb_y<-length(grid_y)
  ntot<-nb_x*nb_y
  syst_p<-matrix(0,nb_x,nb_y)
  dire_p<-matrix(0,nb_x,nb_y)
  syst_a<-matrix(0,nb_x,nb_y)
  dire_a<-matrix(0,nb_x,nb_y)
  bounds<-matrix(0,nb_x,nb_y)
  dir_bounds<-matrix(0,nb_x,nb_y)
  for (k in 1:nb_data)
      {
       m<-matrix(0,nb_x,nb_y)
      n<-matrix(0,nb_x,nb_y)
      o<-matrix(0,nb_x,nb_y)
      p<-matrix(0,nb_x,nb_y)
      .C("systemic_allelic_dom_rg",as.double(m), as.double(n), as.double(o), as.double(p), as.double(coord), as.double(ranks), as.double(genetic),as.double(k), as.integer(k), as.double(grid_x), as.double(grid_y), as.integer(nb_x), as.integer(nb_y), as.integer(nb_ind_mir),as.double(cvx_mat),PACKAGE="wombsoft")->res
      syst_p<-matrix(res[[1]],nb_x,nb_y)
      dire_p<-matrix(res[[2]],nb_x,nb_y)
      syst_a<-matrix(res[[3]],nb_x,nb_y)
      dire_a<-matrix(res[[4]],nb_x,nb_y)
      n<-length(syst_p[syst_p==0])
		  l_quant_p<-quantile(syst_p,probs=c(n/ntot+(1-pB)*(1-n/ntot)))
		  l_quant_a<-quantile(syst_a,probs=c(n/ntot+(1-pB)*(1-n/ntot)))
      for (i in 1:nb_x)
		    {
		    for (j in 1:nb_y)
			   {
			   if (syst_p[i,j]>l_quant_p)
			     {
			     bounds[i,j]<-bounds[i,j]+1
			     dir_bounds[i,j]<-dir_bounds[i,j]+2*dire_p[i,j]
			     }
	       if (syst_a[i,j]>l_quant_a)
			     {
			     bounds[i,j]<-bounds[i,j]+1
			     dir_bounds[i,j]<-dir_bounds[i,j]+2*dire_a[i,j]
			     }
	       }
        }
        rm(m)
        rm(n)
        }
  list(bounds,dir_bounds,pB)
  }
  
 
