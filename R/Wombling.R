####################### DOMINANT ###############################################

DataDominant<-function(input_file,conversion,nb_x,nb_y,output_coords="coord_km.txt")   # data avec nom de l'indiv en 1ère colonne !!!!!!
  {
  input<-read.table(input_file)
  coord<-as.matrix(input[,2:3])
  genetic<-as.matrix(input[,4:ncol(input)])
  if (conversion==1)
    {
    coord<-Conversion(coord)
    } 
  write.table(coord,output_coords)
  min_x<-min(coord[,1])
  max_x<-max(coord[,1])
  min_y<-min(coord[,2])
  max_y<-max(coord[,2])
  grid_x<-seq(min_x,max_x,length=nb_x)
  grid_y<-seq(min_y,max_y,length=nb_y)
  grid<-list(grid_x,grid_y)
  cvx_vertices<-ConvexHull(coord)
  cvx_matrix<-matrix(0,nb_x,nb_y)
  for (i in 1:nb_x)
    {
    for (j in 1:nb_y)
      {
      if (InConvexHull(cvx_vertices,grid_x[i],grid_y[j])==1)
        {cvx_matrix[i,j]<-1}
      }
    }
  list(coord,genetic,grid,cvx_vertices,cvx_matrix,nrow(coord))
  }

MirrorDominant<-function(data,m)    # data = résultat de la fonction Data (donc contient direct plein de trucs)
  {
  grid<-data[[3]]
  cvx<-as.matrix(data[[4]])
  coord<-data[[1]]
  nb_ind<-nrow(coord)
  genetic<-as.matrix(data[[2]])
  plot(coord)
  bord<-Border(grid,cvx)
  CopyDominant(coord,genetic,grid,bord,cvx,m)->res
  list(res[[1]],res[[2]],grid,cvx,data[[5]],nb_ind)               # renvoie le même vecteur, mais avec plein de gens en plus !
  }

WomblingDominant<-function(data,h,output_syst="syst.txt",output_dir="dir.txt")
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
  .C("wombling_aflp",as.double(genetic),as.double(coord),as.integer(nb_indiv_mir),as.integer(nb_data),as.double(grid[[1]]),as.integer(nb_x),as.double(grid[[2]]),as.integer(nb_y),as.double(cvx_mat),as.double(h),as.double(syst),as.double(dire),PACKAGE="wombsoft")->res
  systt<-res[[11]]
  dirr<-res[[12]]  
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
  
  
####################### CO-DOMINANT ############################################

DataCodominant<-function(input_file,conversion,nb_x,nb_y,output_coords="coord_km.txt") ## idem data a une colonne "nom de l'indiv"
  {
  dat<-read.table(input_file)
  coord<-as.matrix(dat[,2:3])
  genetic<-as.matrix(dat[,4:ncol(dat)])
  if (conversion==1)
    {
    coord<-Conversion(coord)
    }
  write.table(coord,output_coords)
  genetic_encoded<-DataEncodingCodominant(genetic) 
  min_x<-min(coord[,1])
  max_x<-max(coord[,1])
  min_y<-min(coord[,2])
  max_y<-max(coord[,2])
  grid_x<-seq(min_x,max_x,length=nb_x)
  grid_y<-seq(min_y,max_y,length=nb_y)
  grid<-list(grid_x,grid_y)
  cvx_vertices<-ConvexHull(coord)
  cvx_matrix<-matrix(0,nb_x,nb_y)
  for (i in 1:nb_x)
    {
    for (j in 1:nb_y)
      {
      if (InConvexHull(cvx_vertices,grid_x[i],grid_y[j])==1)
        {cvx_matrix[i,j]<-1}
      }
    }
  list(coord,genetic_encoded,grid,cvx_vertices,cvx_matrix,nrow(coord))
  }

  
 MirrorCodominant<-function(data,m)    # data = résultat de la fonction Data (donc contient direct plein de trucs)
  {
  grid<-data[[3]]
  cvx<-data[[4]]
  coord<-data[[1]]
  nb_ind<-nrow(coord)
  genetic<-data[[2]]
  plot(coord)
  bord<-Border(grid,cvx)
  CopyCodominant(coord,genetic,grid,bord,cvx,m)->res
  list(res[[1]],res[[2]],grid,cvx,data[[5]],nb_ind)               # renvoie le même vecteur, mais avec plein de gens en plus !
  }
 
WomblingCodominant<-function(data,h,output_syst="syst.txt",output_dir="dir.txt")
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
  .C("wombling_codo",as.double(genetic),as.double(coord),as.integer(nb_ind_mir),as.integer(nb_data),as.integer(nb_all_col),as.double(grid[[1]]),as.integer(nb_x), as.double(grid[[2]]), as.integer(nb_y),as.double(cvx_mat),as.double(h),as.double(syst),as.double(dire),PACKAGE="wombsoft")->res
  systt<-res[[12]]
  dirr<-res[[13]]  
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
  
  
############################## QUANTITATIF #####################################
DataContingency<-function(input_file,conversion,nb_x,nb_y,output_coords="coord_km.txt")   # data avec nom de l'indiv en 1ère colonne !!!!!!
  {
  input<-read.table(input_file)
  coord<-as.matrix(input[,2:3])
  genetic<-as.matrix(input[,4:ncol(input)])
  if (conversion==1)
    {
    coord<-Conversion(coord)
    } 
  write.table(coord,output_coords)
  min_x<-min(coord[,1])
  max_x<-max(coord[,1])
  min_y<-min(coord[,2])
  max_y<-max(coord[,2])
  grid_x<-seq(min_x,max_x,length=nb_x)
  grid_y<-seq(min_y,max_y,length=nb_y)
  grid<-list(grid_x,grid_y)
  cvx_vertices<-ConvexHull(coord)
  cvx_matrix<-matrix(0,nb_x,nb_y)
  for (i in 1:nb_x)
    {
    for (j in 1:nb_y)
      {
      if (InConvexHull(cvx_vertices,grid_x[i],grid_y[j])==1)
        {cvx_matrix[i,j]<-1}
      }
    }
  list(coord,genetic,grid,cvx_vertices,cvx_matrix,nrow(coord))
  }

MirrorContingency<-function(data,m)    # data = résultat de la fonction Data (donc contient direct plein de trucs)
  {
  grid<-data[[3]]
  cvx<-data[[4]]
  coord<-data[[1]]
  nb_ind<-nrow(coord)
  genetic<-data[[2]]
  plot(coord)
  bord<-Border(grid,cvx)
  CopyDominant(coord,genetic,grid,bord,cvx,m)->res
  list(res[[1]],res[[2]],grid,cvx,data[[5]],nb_ind)               # renvoie le même vecteur, mais avec plein de gens en plus !
  }

WomblingContingency<-function(data,h,output_syst="syst.txt",output_dir="dir.txt")
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
  .C("wombling_quantit",as.double(genetic), as.double(coord), as.integer(nb_indiv_mir), as.integer(nb_data),as.double(grid[[1]]),as.integer(nb_x), as.double(grid[[2]]), as.integer(nb_y), as.double(cvx_mat), as.double(h), as.double(syst), as.double(dire),PACKAGE="wombsoft")->res
  systt<-res[[11]]
  dirr<-res[[12]]  
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

  
################################################################################  
################################# UTILITAIRES ##################################  
################################################################################

################################ Enveloppe Convexe #############################  

################################ ConvexHull ####################################

ConvexHull<-function(C)
	{
	pol<-PolarSort(C)
	n<-length(pol)
	i<-1
	while ((i)<=n)
		{
		if (TurnLeft(C[pol[i],1:2],C[pol[i+1-floor(i/n)*n],1:2],C[pol[i+2-floor((i+1)/n)*n],1:2])==1)
			{
			i<-i+1
			}
		else
			{
			pol<-pol[-((i+1)-floor(i/n)*n)]
			if (i!=1)
				{
				i<-i-1
				}
			}
		n<-length(pol)
		}
#	print(C[pol,1:2])
	points(C[pol,1],C[pol,2],pch=19,col='red')
	C[pol,1:2]
	}

################################ InConvexHull #######################################

InConvexHull<-function(cvx,x,y)
	{
	n<-nrow(cvx)
	i<-1
	b<-0
	points(x,y,pch=19,col='lightgrey')
	while (i<=n)
		{
		if (((x>=cvx[i,1])&&(x<=cvx[i+1-floor(i/n)*n,1]))|((x<=cvx[i,1])&&(x>=cvx[i+1-floor(i/n)*n,1])))
			{
			if (TurnLeft(cvx[i,],cvx[i+1-floor(i/n)*n,],c(x,y)))
				{
				if (b==1)
					{
					points(x,y,col='maroon')
					return(1)
					i<-n+1				
					}
				b<-1
				}
			else
				{
				i<-n+1
				return(0)
				}
			}
		i<-i+1
		}
	return(0)
	}

################################ InfPolar ###########################################

InfPolar<-function(A,B)	
	{
	if (A[2]<B[2])
		{
		return(1)
		}
	if (A[2]>B[2])
		{
		return(0)
		}
	if (A[2]==B[2])
		{
		if (A[1]>B[1])
			{
			return(1)
			}
		else
			{
			return(0)
			}
		}
	}

################################ PolarSort ##########################################

PolarSort<-function(C)				# range les points par ordre polaire dont l'origine sera donnée par un point qcq et la demi-droite d'org. 
	{					# par "le" point d'abs max
	
	plot(C[,1],C[,2],pch=19,xlab='x-axis',ylab='y-axis')
	title(main="individuals positions and grid points")
	n<-nrow(C)
	i0<-which.max(C[,1])
	A<-as.double(C[i0,1:2])			# (OA) est la demi-droite d'origine des angles du repère polaire
	x0<-0
	y0<-0
	for (i in 1:4)
		{
		x0<-x0+C[i,1]
		y0<-y0+C[i,2]
		}
	O<-c(x0/4,y0/4)
	points(A[1],A[2],col='blue')
	points(O[1],O[2],col='green')			# origine du repère polaire
	u<-as.double(c((A[1]-O[1])/sqrt((A[1]-O[1])^2+(A[2]-O[2])^2),(A[2]-O[2])/sqrt((A[1]-O[1])^2+(A[2]-O[2])^2)))
	coord_pol<-matrix(0,n,2)
	for (i in 1:n)
		{
		x<-C[i,1]
		y<-C[i,2]
		r<-sqrt((x-O[1])^2+(y-O[2])^2)
		coord_pol[i,1]<-r
		theta<-((x-O[1])*u[1]+(y-O[2])*u[2])/r
		if (theta>1)
			{theta<-1}
		if (theta<(-1))
			{theta<-(-1)}
		coord_pol[i,2]<-as.double(acos(theta))
		if (((y-O[2])*u[1]-(x-O[1])*u[2])<0)
			{
			coord_pol[i,2]<-2*pi-coord_pol[i,2]
			}
		}
	Quicksort(1:n,coord_pol)
	}

################################ Quicksort ############################################

Quicksort<-function(tab,liste_points)		##### ATTENTION : on classe tab par ordre polaire croissant des points de la liste indicée par tab !!!
	{
	n<-length(tab)
	if (n>1)
		{
		tab2<-rep(0,n)
		v<-n+1
		while (v==(n+1))
			{		
			v<-floor(n*runif(1))+1
			}
		min<-1
		max<-n
		for (i in 1:n)
			{
			if (i!=v)
				{
				if (InfPolar(as.double(liste_points[tab[i],]),as.double(liste_points[tab[v],]))==1)
					{
					tab2[min]<-tab[i]
					min<-min+1
					}
				else
					{
					tab2[max]<-tab[i]
					max<-max-1
					}
				}
			}
		if (min==1)
			{
			return(c(tab[v],Quicksort(tab2[(max+1):n],liste_points)))
			}
		else
			{
			if (max==n)
				{
				return(c(Quicksort(tab2[1:(min-1)],liste_points),tab[v]))
				}
			else
				{
				return(c(Quicksort(tab2[1:(min-1)],liste_points),tab[v],Quicksort(tab2[(max+1):n],liste_points)))
				}
			}
		}
	if (n==1)
		{
		return(tab)
		}
	if (n==0)
		{
		return(c())
		}
	}

################################ TurnLeft ###########################################

TurnLeft<-function(p,q,r)
	{
	if ((q[1]-p[1])*(r[2]-p[2])-(q[2]-p[2])*(r[1]-p[1])<=0)
		{return(0)}
	else
		{return(1)}
	}

################################################################################  
################################ DataEncoding ##################################
  
DataEncodingCodominant<-function(D)
	{                                     # ATTENTION : D ne contient pas de coordonnées !
	t<-dim(D)
	nb_indiv<-t[1]
	nb_data<-t[2]
	nb_tot<-0
	DD<-matrix(0,2*nb_indiv,nb_data/2)
	for (i in 1:(nb_indiv))
		{
		for (j in 1:(nb_data/2))
			{
			DD[i,j]<-D[i,2*j-1]
			DD[i+nb_indiv,j]<-D[i,2*j]
			}
		}
	DataEncoding(DD)
	}


DataEncoding<-function(D)			
	{                                   # ATTENTION : D ne contient pas de coordonnées !
	t<-dim(D)
	nb_indiv<-t[1]
	nb_data<-t[2]
	num_alleles<-rep(0,nb_data)
	nb_tot<-0	
	Z<-matrix(0,nb_indiv,nb_data)
	for (j in 1:nb_data)
		{
		alleles<-ListValues(D[,j])
		num_alleles[j]<-length(alleles)
		for (i in 1:nb_indiv)
			{
			if (D[i,j]!=-1)
				{
				b<-0
				k<-1
				while (b==0)
					{
					if (D[i,j]==alleles[k])
						{
						Z[i,j]<-nb_tot+k
						b<-1
						}
					k<-k+1
					}
				}
			else
				{
				Z[i,j]<-0
				}
			}
		nb_tot<-nb_tot+num_alleles[j]
		}
	Z
	}


ListValues<-function(Z)
	{
	l<-length(Z)
	list_tmp<-c(-10)
	for (i in 1:l)
		{
		b<-0
		k<-1
		tmp<-length(list_tmp)
		if (Z[i]!=-1)
			{
			while ((b==0)&(k<=tmp))
				{
				if (Z[i]==list_tmp[k])
					{
					b<-1
					}
				k<-k+1
				}
			if (b==0)
				{
				list_tmp<-c(list_tmp,Z[i])
				}
			}
		}
	tmp<-length(list_tmp)
	list_tmp[2:tmp]
	}


###################### conversion ##############################################
Conversion<-function(D)
	{
	R<-6378			
	pi<-3.14159279
	n<-nrow(D)
	C<-matrix(0,n,2)	
	theta<-0
	phi<-0
	for (i in 1:n)
		{
		theta<-theta+D[i,1]
		phi<-phi+D[i,2]
		}
	theta<-theta/n
	phi<-phi/n
	for (i in 1:n)
		{
		theta_<-D[i,1]
		phi_<-D[i,2]
		C[i,1]<-R*(theta_-theta)*2*pi/360*sin((phi)*2*pi/360)
		C[i,2]<-R*(phi_-phi)*2*pi/360
		}
	C
	}

################### plotDir ####################################################
PlotDir<-function(theta,x,y,length)
	{
	c(x-(length/2)*cos(theta),y-(length/2)*sin(theta),x+(length/2)*cos(theta),y+(length/2)*sin(theta))
	}

################################################################################
################## TOUT MIROIR #################################################

Border<-function(grid,cvx)
	{
	nx<-length(grid[[1]])
	ny<-length(grid[[2]])
	l_bottom<-c()
	points_x<-c()
	l_top<-c()
	l_side<-c()
	points_y<-c()
	for (j in 1:(ny-1))
		{
		if (InConvexHull(cvx,grid[[1]][2],grid[[2]][j])==1)
			{
			l_side<-c(l_side,list(c(grid[[1]][2],grid[[2]][j])))
			points_x<-c(points_x,grid[[1]][2])
			points_y<-c(points_y,grid[[2]][j])
			}
		if (InConvexHull(cvx,grid[[1]][nx-1],grid[[2]][j])==1)
			{
			l_side<-c(l_side,list(c(grid[[1]][nx-1],grid[[2]][j])))
			points_x<-c(points_x,grid[[1]][nx-1])
			points_y<-c(points_y,grid[[2]][j])
			}
		}
	for (i in 3:(nx-1))
		{
		chg<-1
		x0<-grid[[1]][i]
		for (j in 1:(ny-1))
			{
			y0<-grid[[2]][j]
			if (InConvexHull(cvx,x0,y0)==1)
				{
				if (chg==1)
					{
					l_bottom<-c(l_bottom,list(c(x0,y0,BottomSide(cvx,x0,y0))))
					points_x<-c(points_x,grid[[1]][i])
					points_y<-c(points_y,grid[[2]][j])
					chg<-0
					}
				if (InConvexHull(cvx,x0,grid[[2]][j+1])==0)
					{
					l_top<-c(l_top,list(c(x0,y0,TopSide(cvx,x0,y0))))
					points_x<-c(points_x,grid[[1]][i])
					points_y<-c(points_y,grid[[2]][j])
					}
				}
			}
		}
	points(points_x,points_y,col='green')
	list(l_bottom,l_top,l_side)
	}

################################ CopyCodominant ######################################

CopyCodominant<-function(C,Z,grid,border,cvx,h)
	{
	bot<-border[[1]]
	top<-border[[2]]
	side<-border[[3]]
	nb_indiv<-nrow(C)
	nb_vertices<-nrow(cvx)
	CC<-C
	Z1<-as.matrix(Z[1:nb_indiv,])
	Z2<-as.matrix(Z[(nb_indiv+1):(2*nb_indiv),])
	nx<-length(grid[[1]])
	bot_side<-bot[[1]][3]
	ind<-1:(nb_indiv+1)
	ind_bot<-c()
	i<-1
	while (i<length(bot))
		{
		b<-bot[[i]]
		bot_side<-b[3]
		ind_bot<-c(nb_indiv+1)
		while ((b[3]==bot_side)&(i<length(bot)))
			{
			for (k in ind[-ind_bot])
				{
				x<-C[k,1]
				y<-C[k,2]
				if ((b[1]-x)^2+(b[2]-y)^2<=h^2)
					{
					ind_bot<-c(ind_bot,k)
					ww<-Symmetry(c(x,y),cvx[bot_side,1:2],cvx[bot_side+1-floor(bot_side/nb_vertices)*nb_vertices,1:2])
					CC<-rbind(CC,ww)
					Z1<-rbind(Z1,as.vector(Z1[k,]))
					Z2<-rbind(Z2,as.vector(Z2[k,]))
					}
				}
			i<-i+1
			b<-bot[[i]]
			}
		}
	top_side<-top[[1]][3]
	i<-1
	ind_top<-c()
	while (i<length(top))
		{
		b<-top[[i]]
		top_side<-b[3]
		ind_top<-c(nb_indiv+1)
		while ((b[3]==top_side)&(i<length(top)))
			{
			for (k in ind[-ind_top])
				{
				x<-C[k,1]
				y<-C[k,2]
				if ((b[1]-x)^2+(b[2]-y)^2<=h^2)
					{
					ind_top<-c(ind_top,k)
					ww<-Symmetry(c(x,y),cvx[top_side,1:2],cvx[top_side+1-floor(top_side/nb_vertices)*nb_vertices,1:2])
					CC<-rbind(CC,ww)
					Z1<-rbind(Z1,as.vector(Z1[k,]))
					Z2<-rbind(Z2,as.vector(Z2[k,]))
					}
				}
			i<-i+1
			b<-top[[i]]
			}
		}
	for (i in 1:length(side))
		{
		ind_side<-c(nb_indiv+1)
		s<-side[[i]]
		for (k in ind[-ind_side])
			{
			x<-C[k,1]
			y<-C[k,2]
			if ((s[1]-x)^2+(s[2]-y)^2<=h^2)
				{
				ind_side<-c(ind_side,k)
				ww<-Symmetry(c(x,y),c(s[1],0),c(s[1],1))
				CC<-rbind(CC,ww)
				Z1<-rbind(Z1,as.vector(Z1[k,]))
				Z2<-rbind(Z2,as.vector(Z2[k,]))
				}
			}
		}
	points(CC[(nb_indiv+1):nrow(CC),1],CC[(nb_indiv+1):nrow(CC),2])
	list(CC,rbind(Z1,Z2))
	}

################################ CopyDominant ###########################################

CopyDominant<-function(C,Z,grid,border,cvx,h)
	{
	bot<-border[[1]]
	top<-border[[2]]
	side<-border[[3]]
	nb_indiv<-nrow(C)
	nb_vertices<-nrow(cvx)
	CC<-C
	Z<-as.matrix(Z)
	nx<-length(grid[[1]])
	bot_side<-bot[[1]][3]
	ind<-1:(nb_indiv+1)
	ind_bot<-c()
	i<-1
	while (i<length(bot))
		{
		b<-bot[[i]]
		bot_side<-b[3]
		ind_bot<-c(nb_indiv+1)
		while ((b[3]==bot_side)&(i<length(bot)))
			{
			for (k in ind[-ind_bot])
				{
				x<-C[k,1]
				y<-C[k,2]
				if ((b[1]-x)^2+(b[2]-y)^2<=h^2)
					{
					ind_bot<-c(ind_bot,k)
					ww<-Symmetry(c(x,y),cvx[bot_side,1:2],cvx[bot_side+1-floor(bot_side/nb_vertices)*nb_vertices,1:2])
					CC<-rbind(CC,ww)
					Z<-rbind(Z,as.vector(Z[k,]))
					}
				}
			i<-i+1
			b<-bot[[i]]
			}
		}
	top_side<-top[[1]][3]
	i<-1
	ind_top<-c()
	while (i<length(top))
		{
		b<-top[[i]]
		top_side<-b[3]
		ind_top<-c(nb_indiv+1)
		while ((b[3]==top_side)&(i<length(top)))
			{
			for (k in ind[-ind_top])
				{
				x<-C[k,1]
				y<-C[k,2]
				if ((b[1]-x)^2+(b[2]-y)^2<=h^2)
					{
					ind_top<-c(ind_top,k)
					ww<-Symmetry(c(x,y),cvx[top_side,1:2],cvx[top_side+1-floor(top_side/nb_vertices)*nb_vertices,1:2])
					CC<-rbind(CC,ww)
					Z<-rbind(Z,as.vector(Z[k,]))
					}
				}
			i<-i+1
			b<-top[[i]]
			}
		}
	for (i in 1:length(side))
		{
		ind_side<-c(nb_indiv+1)
		s<-side[[i]]
		for (k in ind[-ind_side])
			{
			x<-C[k,1]
			y<-C[k,2]
			if ((s[1]-x)^2+(s[2]-y)^2<=h^2)
				{
				ind_side<-c(ind_side,k)
				ww<-Symmetry(c(x,y),c(s[1],0),c(s[1],1))
				CC<-rbind(CC,ww)
				Z<-rbind(Z,as.vector(Z[k,]))
				}
			}
		}
	points(CC[(nb_indiv+1):nrow(CC),1],CC[(nb_indiv+1):nrow(CC),2])
	list(CC,Z)
	}

################################ TopSide ############################################

TopSide<-function(cvx,x,y)
	{
	n<-nrow(cvx)
	i<-1
	b<-0
	while (i<=n)
		{
		if ((x<cvx[i,1])&&(x>=cvx[i+1-floor(i/n)*n,1]))
			{
			return(i)
			}
		i<-i+1
		}
	}

################################ BottomSide #########################################

BottomSide<-function(cvx,x,y)
	{
	n<-nrow(cvx)
	i<-1
	b<-0
	while (i<=n)
		{
		if ((x>=cvx[i,1])&&(x<cvx[i+1-floor(i/n)*n,1]))
			{
			return(i)
			}
		i<-i+1
		}
	}
	

################################ Symmetry ###########################################

Symmetry<-function(a,b,c)				# symétrique de A par rapport à (BC)
	{
	if (b[1]!=c[1])
		{
		u<-(b[2]-c[2])/(b[1]-c[1])
		v<-(a[1]+(a[2]-b[2]+u*b[1])*u)/(1+u^2)
		xx<-as.double(2*v-a[1])
		yy<-as.double(2*b[2]+2*u*(v-b[1])-a[2])
		}
	else
		{
		yy<-as.double(a[2])
		xx<-as.double(2*b[1]-a[1])
		}
	c(xx,yy)
	}
  
