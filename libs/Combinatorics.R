library(glmnet)
libs_path<-file.path("..","libs")
library(lars)
library(polynom)
library(pls)
source(file.path(libs_path,"Create_synthetic_datasets.R"))
source(file.path(libs_path,"helper_functions.R"))


assert <- function(condition, message) {
  if (!condition) stop(message)
}


get_range4<- function(x,l1=21,l2=14,l3=2,l4=3) #returns the inddices from the range of x
{ assert(x<=l1+l2+l3+l4, "x should be in correct range")
  if (x<=l1)
{return(c(1:l1))}
  if (x<=l1+l2)
  {return(c( (l1+1) : (l1+l2) ))}
     if (x<=l1+l2+l3)
     {return(c( (l1+l2+1) : (l1+l2+l3) )) }
        return(c( (l1+l2+l3+1):(l1+l2+l3+l4) ) )
}

#get_range4(23)

get_ranges4<-function(l1=21,l2=14,l3=2,l4=3)
{ l_main<-l1+l2+l3+l4
  l_theta<- l1*(l2+l3+l4)+l2*(l3+l4) +l3*l4
  l_psi<- l1*l2*(l3+l4) +l3*l4*(l1+l2)
  l_phi<-l1*l2*l3*l4
  
range_main<-c(1: l_main )
range_theta<-c( (l_main+1) : (l_main+l_theta) )
range_psi<-c( (l_main + l_theta +1) : (l_main + l_theta + l_psi) )
range_phi<-c( (l_main + l_theta + l_psi+1): (l_main + l_theta + l_psi + l_phi))

return(list(range_main, range_theta, range_psi, range_phi))}

print(get_ranges4(l1=2,l2=1,l3=2,l4=3))


#################### FUNCTIONS 2-way combinatorics ################################

#Function that transform matrix theta_hat in vec_theta 
get_theta_vec_2way4<-function(Theta_hat, l1=21,l2=14,l3=2,l4=3)
  
{range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))
range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
counter<-0
vec_theta<-array(0,l1*(l2+l3+l4)+l2*(l3+l4) +l3*l4)

## case 1 :a with b or c or d}
for (i in range1)
{for (j in c(range2,range3, range4))
{counter<-counter+1
vec_theta[counter]<-(Theta_hat[i,j] +Theta_hat[j,i])/2
}}

## case 2: b with c or d
for (i in range2)
{for (j in c(range3, range4) )
{counter<-counter+1
vec_theta[counter]<-(Theta_hat[i,j] +Theta_hat[j,i])/2
}}


## case 3: c with d
for (i in range3)
{for (j in  range4 )
{counter<-counter+1
vec_theta[counter]<-(Theta_hat[i,j] +Theta_hat[j,i])/2
}}

return(vec_theta)
}

#get_theta_vec_2way4(matrix(c(c(1,2,3,4,5,6),c(1,2,3,4,5,6),c(1,2,3,4,5,6),c(1,2,3,4,5,6),c(1,2,3,4,5,6),c(1,2,3,4,5,6)), nrow=6), l1=1,l2=1,l3=1,l4=3)



get_theta_from_theta_vec_2way4<-function(vec_theta,l1=21,l2=14, l3=2, l4=3) #get theta matrix from theta vec
{ counter<-1
range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))
range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
Theta_hat<-matrix(0,nrow=l1+l2+l3+l4, ncol=l1+l2+l3+l4)

## case 1 :a with b or c or d
for (i in range1)
{for (j in c(range2,range3, range4))
{
  Theta_hat[i,j]<-vec_theta[counter]
  Theta_hat[j,i]<-vec_theta[counter] 
  counter<-counter+1
}}

## case 2: b with c or d
for (i in range2)
{for (j in c(range3, range4) )
{
  Theta_hat[i,j]<-vec_theta[counter] 
  Theta_hat[j,i]<-vec_theta[counter] 
  counter<-counter+1
}}

## case 2: c with d
for (i in range3)
{for (j in  range4 )
{
  Theta_hat[i,j]<-vec_theta[counter] 
  Theta_hat[j,i]<-vec_theta[counter] 
  counter<-counter+1
}}

assert(counter==l1*l2+l2*l3+l3*l1 +l4*(l1+l2+l3)+1, 'smth wrong with counter')
return(Theta_hat)
}


get_theta_from_theta_vec_2way4(vec_theta = c(1,1,1,1,10,10,10,10,2,2,2,3,3),l1=2,l2=1, l3=1, l4=2)



##position in vector form to position in matrix form 2way
get_position_vec_from_theta_matrix4<- function(position_tuple, l1=21,l2=14,l3=2,l4=3) ## takes into account / works only for possible combinations!!!!
{ x<-position_tuple[1]
y<-position_tuple[2]

range_x<-get_range4(x,l1=l1,l2=l2,l3=l3, l4=l4)
range_y<-get_range4(y,l1=l1,l2=l2,l3=l3,l4=l4)


assert(x<= l1+l2+l3+l4, "x should be <=l1+l2")
assert(x<y, "x<y")
assert(y>l1, 'y should be >l1')


if( all(range_x == c(1:l1)) ==TRUE ) #ab or ac or ad
{ 
  position_vector<- (x-1)*(l2+l3+l4) +(y-l1)  }


if( all ( range_x == c( (l1+1): (l1+l2) ) ) == TRUE )  #bc or bd

{position_vector<-l1*(l2+l3+l4) + (x-l1-1)*(l3+l4) + y- (l1+l2)  } 


if( all ( range_x == c( (l1+l2+1): (l1+l2+l3) ) ) == TRUE )  #cd
{position_vector<-l1*(l2+l3+l4) + l2*(l3+l4) + (x-l1-l2-1)*l4 + y - (l1+l2+l3)  } 
return(position_vector)
}

get_position_vec_from_theta_matrix4(c(22,36))


get_beta_vec_2way4<-function(beta,l1,l2,l3,l4, gamma, only_beta = FALSE) #beta is beta_main
{range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))
range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
beta_vec2way<- array(0, dim = l1*(l2+l3+l4) +l2*l3 +l2*l4 +l3*l4 )
counter<-1

for (i in range1){ #ab ac ad
  for (j in  c(range2, range3, range4) ){
    beta_vec2way[counter]<-beta[i]*beta[j]
    counter<-counter+1}}

for (i in range2){ #bc bd
  for (j in  c( range3, range4 ) ){
    beta_vec2way[counter]<-beta[i]*beta[j]
    counter<-counter+1}}

for (i in range3){ #cd
  for (j in   range4  ){
    beta_vec2way[counter]<-beta[i]*beta[j]
    counter<-counter+1}}


assert(counter==l1*(l2+l3+l4)+l2*l3 +l2*l4 + l3*l4 +1)
if (only_beta==FALSE)
{beta_vec2way<-beta_vec2way*gamma}
return(beta_vec2way)
}

get_beta_vec_2way4(beta=c(1,2,3,4, 5), gamma=1, l1=1, l2=1, l3=1, l4=2)


#gets ls of positions in vector form  from ls positions in matrix
get_positions_2way4<-function(ls_positions,l1=21,l2=14,l3=2,l4=3){
  all_positions<-c()
  for (tuple in ls_positions)
  {pos<-get_position_vec_from_theta_matrix4(position_tuple = tuple, l1=l1, l2=l2, l3=l3, l4=l4)
  all_positions<-c(all_positions, pos)}
  return(all_positions)}

x<-get_positions_2way4(ls_positions=list(c(1,22), c(36,38)),l1=21,l2=14,l3=2,l4=3)



#####################  FUNCTIONS 3-way combinatorics  ####################################


#value of psi_vec at tablepsi_[i,j,k]
psi_value_from_table_position<-function (table,i,j,k)
{return( (table[i,j,k] + table[i,k,j] + table [j,i,k] +table[j,k,i] + table[k,i,j] + table[k,j,i] )/6)}

#position in 3 dim table to vector index: works only for possible combinations
psi_table_position_to_vector_index4<- function(position_tuple,l1=21,l2=14,l3=2,l4=3) ## takes into account / works only for possible combinations!!!!
{
  
 
  
  x<-position_tuple[1]
  y<-position_tuple[2]
  z<-position_tuple[3]
  #print(y)
  
  range_x<-get_range4(x,l1=l1,l2=l2,l3=l3, l4=l4)
  range_y<-get_range4(y,l1=l1,l2=l2,l3=l3,l4=l4)
  range_z<-get_range4(z,l1=l1,l2=l2,l3=l3,l4=l4)
  #print(range_x)
  #print(range_y)
  #print(range_z)
  
  assert(x<= l1+l2+l3+l4, "x should be <=l1+l2")
  assert(x<y, "x<y")
  assert(y>l1, 'y should be >l1')
  assert(y<z)
  
  #print("ok")
  if( all(range_x == c(1:l1)) ==TRUE ) #ab c/d, acd
    #print("ok x")
  { if (all(range_y == c( (l1+1): (l1+l2) )) ==TRUE  ) { #ab c,d
    #print("ok y")
    position_vector<- (x-1)* (l2* (l3+l4) +l3*l4) + (y-l1-1)*(l3+l4) + (z-l1-l2)}
    if (all(range_y == c( (l1+l2+1): (l1+l2+l3) )) ==TRUE  ){  #acd  
      #print("ok acd")
      position_vector<-    (x-1)* (l2* (l3+l4) +l3*l4) + l2*(l3+l4) + (y-l1-l2-1)*l4 + (z-l1-l2-l3)}
    
    
    }
  #print("ok2")
  
  if( all ( range_x == c( (l1+1): (l1+l2) ) ) == TRUE )  #bcd
    
  {position_vector<-l1*(l2*(l3+l4) +l3*l4) + (x-l1-1)*(l3*l4) +( y- (l1+l2+1) ) *l4 +  z-(l1+l2+l3)  } 
  
  return(position_vector)

  }


#psi_table_position_to_vector_index4(position_tuple=c(1, 23,40),l1=21,l2=14,l3=2,l4=3)


#get psi vec from 3dim table
get_psi_vec4<-function(psi, l1=21,l2=14,l3=2,l4=3)
{
  assert(all(dim(psi)==l1+l2+l3+l4), "Dimensions are not ok")
  
  psi_vec<-array(0, dim=c(l1*l2*l3+l1*l2*l4+l1*l3*l4+l2*l3*l4) )
  range1<-c(1:l1)
  range2<-c((l1+1):(l1+l2))
  range3<-c((l1+l2+1):(l1+l2+l3))
  range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
  
  #print("start for")
  
  for (i in range1) { #abc
    for (j in range2 ) { 
      for (k in range3 ) {  
        psi_vec[psi_table_position_to_vector_index4(c(i,j,k),l1=l1,l2=l2,l3=l3, l4=l4)]<-psi_value_from_table_position(psi,i,j,k)
      }}}
  
  for (i in range1) { #abd
    for (j in range2 ) { 
      for (l in range4 ) {  
        psi_vec[psi_table_position_to_vector_index4(c(i,j,l),l1=l1,l2=l2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,i,j,l)
      }}}
  
  for (i in range1) { #abd
    for (k in range3 ) { 
      for (l in range4 ) {  
        psi_vec[psi_table_position_to_vector_index4(c(i,k,l),l1=l1,l2=l2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,i,k,l)
      }}}
  
  for (j in range2) { #abd
    for (k in range3 ) { 
      for (l in range4 ) {  
        psi_vec[psi_table_position_to_vector_index4(c(j,k,l),l1=l1,l2=l2,l3=l3,l4=l4)]<-psi_value_from_table_position(psi,j,k,l)
      }}}
  
  
  
  return(psi_vec)
}

psi=array(1, dim = c(4, 4, 4))
dim(psi)
#get_psi_vec4(psi, l1=1, l2=1, l3=1 , l4=1)

#get psi table from psi_vec 
get_psi_from_psi_vec4<-function(psi_vec,l1=21,l2=14,l3=2,l4=3)
{ range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))
range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
  counter<-1
psi<-array(0,dim=(c(l1+l2+l3+l4,l1+l2+l3+l4,l1+l2+l3+l4)))

for (i in range1) { 
  for (j in range2 ) { 
    for (k in c(range3,range4) ) { 
      #cat("i:", i, ", j:", j, ", k:", k ,'\n')
      psi[i,j,k]<-psi_vec[counter]
      psi[i,k,j]<-psi_vec[counter]
      psi[k,i,j]<-psi_vec[counter]
      psi[k,j,i]<-psi_vec[counter]
      psi[j,i,k]<-psi_vec[counter]
      psi[j,k,i]<-psi_vec[counter]
      counter<-counter+1
    }}} 

for (i in range1) { 
  for (j in range3 ) { 
    for (k in range4 ) { 
      #cat("i:", i, ", j:", j, ", k:", k ,'\n')
      psi[i,j,k]<-psi_vec[counter]
      psi[i,k,j]<-psi_vec[counter]
      psi[k,i,j]<-psi_vec[counter]
      psi[k,j,i]<-psi_vec[counter]
      psi[j,i,k]<-psi_vec[counter]
      psi[j,k,i]<-psi_vec[counter]
      counter<-counter+1
    }}} 

for (i in range2) { 
  for (j in range3 ) { 
    for (k in range4 ) { 
      print(counter)
      #cat("i:", i, ", j:", j, ", k:", k ,'\n')
      psi[i,j,k]<-psi_vec[counter]
      psi[i,k,j]<-psi_vec[counter]
      psi[k,i,j]<-psi_vec[counter]
      psi[k,j,i]<-psi_vec[counter]
      psi[j,i,k]<-psi_vec[counter]
      psi[j,k,i]<-psi_vec[counter]
      counter<-counter+1
    }}} 

assert(counter==1+ l1*l2*(l3+l4) +l3*l4*(l1+l2))
return(psi)
}

#psi<-get_psi_from_psi_vec4(psi_vec=6*c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),l1=1,l2=2,l3=2,l4=2)
#psi_vec<-get_psi_vec4(psi, l1=1,l2=2,l3=2,l4=2)  


get_beta_vec_3way4<-function(beta_2way, delta, l1=21,l2=14,l3=2,l4=3, only_beta=FALSE) ## # beta_2way should be final beta_2way; #only_beta means product of beta_2ways (gamma included) without delta
{beta_vec3way<- array(0, dim = l1*l2*l3 )
counter<-1
range1<-c(1:l1)
range2<-c ( (l1+1): (l1+l2) )
range3<- c( (l1+l2+1): (l1+l2+l3) )
range4 <- c ( (l1+l2+l3+1): (l1+l2+l3+l4) )

#Iterate over possible positions
for (i in range1){ #ab c/d
  for (j in range2 ){
    for (k in c(range3, range4) ){
      
        beta_vec3way[counter]<-beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(i,j),l1=l1,l2=l2,l3=l3)]*
        beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(i,k),l1=l1,l2=l2,l3=l3, l4=l4)]*
        beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(j,k),l1=l1,l2=l2,l3=l3, l4=l4)]
      
      counter<-counter+1}}}

for (i in range1){ #acd 
  for (j in range3 ){
    for (k in  range4 ){
      
      beta_vec3way[counter]<-beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(i,j),l1=l1,l2=l2,l3=l3)]*
        beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(i,k),l1=l1,l2=l2,l3=l3, l4=l4)]*
        beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(j,k),l1=l1,l2=l2,l3=l3, l4=l4)]
      
      counter<-counter+1}}}


for (i in range2){ #bcd 
  for (j in range3 ){
    for (k in  range4 ){
      
        beta_vec3way[counter]<-beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(i,j),l1=l1,l2=l2,l3=l3, l4=l4)]*
        beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(i,k),l1=l1,l2=l2,l3=l3, l4=l4)]*
        beta_2way[get_position_vec_from_theta_matrix4(position_tuple = c(j,k),l1=l1,l2=l2,l3=l3, l4=l4)]
      
      counter<-counter+1}}}


if (only_beta == FALSE)
{beta_vec3way<-beta_vec3way*delta}


assert(counter==l1*l2*(l3+l4) +l3*l4*(l1+l2) +1)

return(beta_vec3way)

}

#get_beta_vec_3way4(beta_2way=c(1,1,1,1,1,1,1,1,2,2,2,3,3),l1=2,l2=1,l3=1,l4=2, delta=1, only_beta=FALSE)



get_positions_3way4<-function(ls_positions, l1=21,l2=14,l3=2,l4=3)
{all_positions<-c()
for (tuple in ls_positions)
{pos<-psi_table_position_to_vector_index4(position_tuple = tuple, l1=l1, l2=l2, l3=l3,l4=l4)
#print(pos)
all_positions<-c(all_positions, pos)}
return(all_positions)}

ls_pos<-list(c(1,3,6), c(1,3,5), c(1,4,7), c(2,4,7), c(3,6,7))
get_positions_3way4(ls_pos, l1=1, l2=2, l3=3,l4=1)



#####################  FUNCTIONS 4-way combinatorics  ####################################

phi_value_from_table_position<-function (table,i,j,k,l)
{return( (table[i,j,k,l] + table[i,j,l,k] + table [i,k,j,l] +table[i,k,l,j] + table[i,l,j,k] + table[i,l,k,j] 
          + table[j,i,k,l] + table[j,i,l,k] + table [j,k,i,l] +table[j,k,l,i] + table[j,l,i,k] + table[j,l,k,i] 
          + table[k,i,j,l] + table[k,i,l,j] + table [k,j,i,l] +table[k,j,l,i] + table[k,l,i,j] + table[k,l,j,i] 
          + table[l,i,j,k] + table[l,i,k,j] + table [l,j,i,k] +table[l,j,k,i] + table[l,k,i,j] + table[l,k,j,i]    
)/24
         )}



#Give position is phi_vec from position in 4d table for phi
phi_table_position_to_vector_index4<- function(position_tuple, l1=21,l2=14,l3=2,l4=3) ## takes into account / works only for possible combinations!!!!
{
  
  x<-position_tuple[1]
  y<-position_tuple[2]
  z<-position_tuple[3]
  t<-position_tuple[4]
  
  assert(x<=l1, "x should be <=l1")
  assert(y<=l1+l2, "y should be <=l1+l2")
  assert(z<= l1+l2+l3, "z should be <= than l1+l2+l3")
  assert(t<=l1+l2+l3+l4, "t should be <= l1+l2+l3+l4")
  
  assert(x>=1, "x should be >=1")
  assert(y>l1, "y should be >l1")
  assert(z>l1+l2, 'z should be >l1+l2')
  assert(t>l1+l2+l3)
  
  
  position_phi<-(x-1)*l2*l3*l4 + (y-l1-1)*l3*l4 + (z-l1-l2-1)*l4 + (t-l1-l2-l3)
  
  return(position_phi)
  
}


#phi_table_position_to_vector_index4(position_tuple = c(1,3,6,10),1,2,3,4)

get_phi_vec4<-function(phi, l1=21,l2=14,l3=2,l4=3)
{
  assert(all(dim(phi)==l1+l2+l3+l4), "Dimensions are not ok")
  
  phi_vec<-array(0, dim=l1*l2*l3*l4 )
  range1<-c(1:l1)
  range2<-c((l1+1):(l1+l2))
  range3<-c((l1+l2+1):(l1+l2+l3))
  range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
  
  #print("start for")
  
  for (i in range1) { #abc
    for (j in range2 ) { 
      for (k in range3 ) {  
        for(l in range4)
        phi_vec[phi_table_position_to_vector_index4(c(i,j,k,l),l1=l1,l2=l2,l3=l3, l4=l4)]<-phi_value_from_table_position(phi,i,j,k,l)
      }}}
  
  return(phi_vec)
}

phi=array(2, dim = c(14, 14, 14, 14))
dim(phi)
#get_phi_vec4(phi, l1=2, l2=3, l3=4 , l4=5) 



#get psi table from psi_vec 
get_phi_from_phi_vec4<-function(phi_vec,l1=21,l2=14,l3=2,l4=3)
{ range1<-c(1:l1)
range2<-c((l1+1):(l1+l2))
range3<-c((l1+l2+1):(l1+l2+l3))
range4<-c((l1+l2+l3+1):(l1+l2+l3+l4))
counter<-1
phi<-array(0,dim=(c(l1+l2+l3+l4,l1+l2+l3+l4,l1+l2+l3+l4, l1+l2+l3+l4)))

for (i in range1) { 
  for (j in range2 ) { 
    for (k in range3) {
      for(l in range4){
      #cat("i:", i, ", j:", j, ", k:", k ,'\n')
      phi[i,j,k,l]<-phi_vec[counter]
      phi[i,j,l,k]<-phi_vec[counter]
      phi[i,k,j,l]<-phi_vec[counter]
      phi[i,k,l,j]<-phi_vec[counter]
      phi[i,l,j,k]<-phi_vec[counter]
      phi[i,l,k,j]<-phi_vec[counter]
      
      phi[j,i,k,l]<-phi_vec[counter]
      phi[j,i,l,k]<-phi_vec[counter]
      phi[j,k,i,l]<-phi_vec[counter]
      phi[j,k,l,i]<-phi_vec[counter]
      phi[j,l,i,k]<-phi_vec[counter]
      phi[j,l,k,i]<-phi_vec[counter]
      
      phi[k,i,j,l]<-phi_vec[counter]
      phi[k,i,l,j]<-phi_vec[counter]
      phi[k,j,i,l]<-phi_vec[counter]
      phi[k,j,l,i]<-phi_vec[counter]
      phi[k,l,i,j]<-phi_vec[counter]
      phi[k,l,j,i]<-phi_vec[counter]
      
      phi[l,i,j,k]<-phi_vec[counter]
      phi[l,i,k,j]<-phi_vec[counter]
      phi[l,j,i,k]<-phi_vec[counter]
      phi[l,j,k,i]<-phi_vec[counter]
      phi[l,k,i,j]<-phi_vec[counter]
      phi[l,k,j,i]<-phi_vec[counter]
      print(phi_vec[counter])
      
      counter<-counter+1
    }}} }
assert(counter==1+ l1*l2*l3*l4,"counter error in phi " )
return(phi)
}

#x<-get_phi_from_phi_vec4(phi_vec=c(1:24),l1=3,l2=2,l3=4,l4=1)
#x_vec<-get_phi_vec4(x,l1=3,l2=2,l3=4,l4=1)




get_beta_vec_4way4<-function(beta_3way, tau, l1=21,l2=14,l3=2,l4=3, only_beta=FALSE) ## # beta_2way should be final beta_2way; #only_beta means product of beta_2ways (gamma included) without delta
{ assert(length(beta_3way)==l1*l2*(l3+l4) + l3*l4*(l1+l2))
  beta_vec4way<- array(0, dim = l1*l2*l3*l4 )
counter<-1
range1<-c(1:l1)
range2<-c ( (l1+1): (l1+l2) )
range3<- c( (l1+l2+1): (l1+l2+l3) )
range4 <- c ( (l1+l2+l3+1): (l1+l2+l3+l4) )

#Iterate over possible positions
for (i in range1){ #ab c/d
  for (j in range2 ){
    for (k in range3 ){
      {for (l in range4){
      
      beta_vec4way[counter]<-beta_3way[psi_table_position_to_vector_index4(position_tuple = c(j,k,l),l1=l1,l2=l2,l3=l3,l4=l4)]*
        beta_3way[psi_table_position_to_vector_index4(position_tuple = c(i,j,k),l1=l1,l2=l2,l3=l3, l4=l4)]*
        beta_3way[psi_table_position_to_vector_index4(position_tuple = c(i,j,l),l1=l1,l2=l2,l3=l3, l4=l4)]*
        beta_3way[psi_table_position_to_vector_index4(position_tuple = c(i,k,l),l1=l1,l2=l2,l3=l3, l4=l4)]
        
      counter<-counter+1}}}}}


if (only_beta == FALSE)
{beta_vec4way<-beta_vec4way*tau}


assert(counter==l1*l2*l3*l4 +1)

return(beta_vec4way)

}

#get_beta_vec_4way4(beta_3way=c(1,2,3,4,5,6,7),l1=1,l2=1,l3=1,l4=2, tau=c(1,2), only_beta=FALSE)


get_positions_4way4<-function(ls_positions, l1=21,l2=14,l3=2,l4=3)
{all_positions<-c()
for (tuple in ls_positions)
{pos<-phi_table_position_to_vector_index4(position_tuple = tuple, l1=l1, l2=l2, l3=l3,l4=l4)
#print(pos)
all_positions<-c(all_positions, pos)}
return(all_positions)}

#ls_pos<-list(c(1,2,4,7),c(1,2,5,7), c(1,3,5,7),  c(1,3,6,7))
#get_positions_4way4(ls_pos, l1=1, l2=2, l3=3,l4=1)




######################################################################################################################################################





















































