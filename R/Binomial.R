# Data c'est list(coord,genetic_encoded,grid,cvx_vertices,cvx_matrix,nrow(coord))

################################################################################
###################### CODOMINANT ##############################################

 
 CandidateBoundariesCodominant<-function(data,h,pB=0.3)
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
      .C("systemic_allelic_codom",as.double(m),as.double(n),as.double(coord),as.double(genetic),as.double(h),as.integer(k),as.integer(ou_alleles[k]),as.double(grid_x),as.double(grid_y),as.integer(nb_x),as.integer(nb_y),as.integer(nb_ind_mir),as.double(cvx_mat),PACKAGE="wombsoft")->res
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
  
  BinomialTestCodominant<-function(data,cbe,pvalue=0.05,output_bounds="bounds.txt",output_dir="direction.txt")
    {
    seuil<-cbe[[3]]
    coord<-as.matrix(data[[1]])
    genetic<-data[[2]]
    grid<-data[[3]]
    grid_x<-grid[[1]]
    grid_y<-grid[[2]]
    cvx_mat<-data[[5]]
    nb_data<-2*ncol(genetic)
    nb_x<-length(grid_x)
    nb_y<-length(grid_y)
    barr<-matrix(-1,nb_x,nb_y)
    dir_bar<-matrix(0,nb_x,nb_y)
    nb_alleles<-max(genetic)
    bounds<-cbe[[1]]
    dir<-cbe[[2]]
    seuil_bin<-qbinom(1-pvalue,nb_alleles,seuil)
    for (i in 1:nb_x)
		 {
		 for (j in 1:nb_y)
		  {
		  if (cvx_mat[i,j]==1)
		    {
		    barr[i,j]<-0
		    if (bounds[i,j]>seuil_bin)
		      {
		      barr[i,j]<-1;
		      dir_bar[i,j]<-dir[i,j]/(2*bounds[i,j])
		      }
        }
      }
      }
    write.table(barr,output_bounds)
    write.table(dir_bar,output_dir)
    
    area<-sqrt((grid[[1]][nb_x]-grid[[1]][1])^2+(grid[[2]][nb_y]-grid[[2]][1])^2)
	 list_seg<-list(c(),c(),c(),c())
	 for (i in 1:((nb_x)/3))
		{
		for (j in 1:((nb_y)/3))
			{
			if (barr[3*i,3*j]==1)
				{
				gr<-PlotDir(as.double(dir_bar[3*i,3*j]),as.double(grid_x[3*i]),as.double(grid_y[3*j]),(4*area/(nb_x)))
				for (k in 1:4)
					{
					list_seg[[k]]<-c(list_seg[[k]],gr[k])
					}
				}
			}
		}
	 image(grid_x,grid_y,barr,zlim=c(0,1),col=c('darkgreen','grey'),xlab='x_coordinates',ylab='y_coordinates')
   title(main='Boundaries')	
	 segments(list_seg[[1]],list_seg[[2]],list_seg[[3]],list_seg[[4]],col='black')
	 points(coord[1:data[[6]],])
   }
 
 
################################################################################ 
######################### DOMINANT et QUANTITATIF ##############################
################################################################################ 
 
 CandidateBoundariesDominant<-function(data,h,pB=0.3)
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
      .C("systemic_allelic_dom",as.double(m), as.double(n), as.double(o), as.double(p), as.double(coord), as.double(genetic),as.double(h), as.integer(k), as.double(grid_x), as.double(grid_y), as.integer(nb_x), as.integer(nb_y), as.integer(nb_ind_mir),as.double(cvx_mat),PACKAGE="wombsoft")->res
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
  
  CandidateBoundariesContingency<-function(data,h,pB=0.3)
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
  syst<-matrix(0,nb_x,nb_y)
  dire<-matrix(0,nb_x,nb_y)
  bounds<-matrix(0,nb_x,nb_y)
  dir_bounds<-matrix(0,nb_x,nb_y)
  for (k in 1:nb_data)
      {
      m<-matrix(0,nb_x,nb_y)
      n<-matrix(0,nb_x,nb_y)
      .C("systemic_allelic_quantit",as.double(m), as.double(n),as.double(coord),as.double(genetic),as.double(h), as.integer(k), as.double(grid_x), as.double(grid_y), as.integer(nb_x), as.integer(nb_y), as.integer(nb_ind_mir),as.double(cvx_mat),PACKAGE="wombsoft")->res
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
 
  
  BinomialTestDominant<-function(data,cbe,pvalue=0.05,output_bounds="bounds.txt",output_dir="direction.txt")
    {
    seuil=cbe[[3]]
    coord<-as.matrix(data[[1]])
    genetic<-data[[2]]
    grid<-data[[3]]
    grid_x<-grid[[1]]
    grid_y<-grid[[2]]
    cvx_mat<-data[[5]]
    nb_x<-length(grid_x)
    nb_y<-length(grid_y)
    barr<-matrix(-1,nb_x,nb_y)
    dir_bar<-matrix(0,nb_x,nb_y)
    nb_alleles<-2*ncol(genetic)
    bounds<-cbe[[1]]
    dir<-cbe[[2]]
    seuil_bin<-qbinom(1-pvalue,nb_alleles,seuil)
    for (i in 1:nb_x)
		 {
		 for (j in 1:nb_y)
		  {
		  if (cvx_mat[i,j]==1)
		    {
		    barr[i,j]<-0
		    if (bounds[i,j]>seuil_bin)
		      {
		      barr[i,j]<-1;
		      dir_bar[i,j]<-dir[i,j]/(2*bounds[i,j])
		      }
        }
      }
      }
    write.table(barr,output_bounds)
    write.table(dir_bar,output_dir)
   
    area<-sqrt((grid[[1]][nb_x]-grid[[1]][1])^2+(grid[[2]][nb_y]-grid[[2]][1])^2)
	 list_seg<-list(c(),c(),c(),c())
	 for (i in 1:((nb_x)/3))
		{
		for (j in 1:((nb_y)/3))
			{
			if (barr[3*i,3*j]==1)
				{
				gr<-PlotDir(as.double(dir_bar[3*i,3*j]),as.double(grid_x[3*i]),as.double(grid_y[3*j]),(4*area/(nb_x)))
				for (k in 1:4)
					{
					list_seg[[k]]<-c(list_seg[[k]],gr[k])
					}
				}
			}
		}
	 image(grid_x,grid_y,barr,zlim=c(0,1),col=c('darkgreen','grey'),xlab='x_coordinates',ylab='y_coordinates')	
	 segments(list_seg[[1]],list_seg[[2]],list_seg[[3]],list_seg[[4]],col='black')
	 title(main='Boundaries')
	 points(coord[1:data[[6]],])
   }
   
 BinomialTestContingency<-function(data,cbe,pvalue=0.05,output_bounds="bounds.txt",output_dir="direction.txt")
    {
    seuil<-cbe[[3]]
    coord<-as.matrix(data[[1]])
    genetic<-data[[2]]
    grid<-data[[3]]
    grid_x<-grid[[1]]
    grid_y<-grid[[2]]
    cvx_mat<-data[[5]]
    nb_x<-length(grid_x)
    nb_y<-length(grid_y)
    barr<-matrix(-1,nb_x,nb_y)
    dir_bar<-matrix(0,nb_x,nb_y)
    nb_alleles<-ncol(genetic)
    bounds<-cbe[[1]]
    dir<-cbe[[2]]
    seuil_bin<-qbinom(1-pvalue,nb_alleles,seuil)
    for (i in 1:nb_x)
		 {
		 for (j in 1:nb_y)
		  {
		  if (cvx_mat[i,j]==1)
		    {
		    barr[i,j]<-0
		    if (bounds[i,j]>seuil_bin)
		      {
		      barr[i,j]<-1;
		      dir_bar[i,j]<-dir[i,j]/(2*bounds[i,j])
		      }
        }
      }
      }
    
    write.table(barr,output_bounds)
    write.table(dir_bar,output_dir)
    area<-sqrt((grid[[1]][nb_x]-grid[[1]][1])^2+(grid[[2]][nb_y]-grid[[2]][1])^2)
	 list_seg<-list(c(),c(),c(),c())
	 for (i in 1:((nb_x)/3))
		{
		for (j in 1:((nb_y)/3))
			{
			if (barr[3*i,3*j]==1)
				{
				gr<-PlotDir(as.double(dir_bar[3*i,3*j]),as.double(grid_x[3*i]),as.double(grid_y[3*j]),(4*area/(nb_x)))
				for (k in 1:4)
					{
					list_seg[[k]]<-c(list_seg[[k]],gr[k])
					}
				}
			}
		}
	 image(grid_x,grid_y,barr,zlim=c(0,1),col=c('darkgreen','grey'))	
	 segments(list_seg[[1]],list_seg[[2]],list_seg[[3]],list_seg[[4]],col='black',xlab='x_coordinates',ylab='y_coordinates')
	 title(main="Boundaries")
	 points(coord[1:data[[6]],])
   }
 
 
