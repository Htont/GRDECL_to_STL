# GRDECL reader

gr_path = "" #Path to directory with file
gr_name = "ECLIPSE100_GRID.GRDECL" #Name of file with Grid
gr_prop = "ECLIPSE100_GRID_PROPS.GRDECL" #Name of file with NTG

tmp_file <- file(paste0(gr_path,gr_name), open = "r")
s <- readLines(tmp_file)
close(tmp_file)

nstr = length(s)

tmp_i <- 0
tmp_i2 <- 0
tmp_i3 <- 0
tmp_i4 <- 0


##### SPECGRID #####

for(i in 1:nstr){
  if(s[i] == "SPECGRID"){
    tmp_row <- strsplit(s[i+1], split = "\\s+")
    
    nx = as.numeric(tmp_row[[1]][1])
    ny = as.numeric(tmp_row[[1]][2])
    nz = as.numeric(tmp_row[[1]][3])
    
    m_grd <- array(dim = c(nx,ny,nz))
    l_grd <- array(dim = c(nx,ny,nz))
    pills <- array(dim = c(nx+1,ny+1,2,3))
    z_coords <- array(dim = c(2*nx,2*ny,2*nz))
    
    tmp_i <- i + 1
    break()
  }
}


##### COORD #####

i_x <- i_y <- 1
flag1 <- F

for(i in tmp_i:nstr){
  
  if(s[i] == "COORD"){
    flag1 <- T
    i <- i + 1
  }
  else if(s[i] == "/" & flag1){
    flag1 <- F
    tmp_i2 <- i + 1
    break()
  }
  else if(flag1){
    tmp_row <- as.numeric(strsplit(s[i], split = "\\s+")[[1]])
    
   if(length(tmp_row) == 6){
      pills[i_x,i_y,1,1:3] <- tmp_row[1:3]
      pills[i_x,i_y,2,1:3] <- tmp_row[4:6]
      i_x <- i_x + 1
      if(i_x > nx + 1){
        i_x <- 1
        i_y <- i_y + 1
        if(i_y > ny + 1){
          print("OK")
        }
      }
    }
  }
}

max_depth <- max(pills[,,2,3])

##### ZCORN #####

i_x <- i_y <- i_z <- 1
flag1 <- F

for(i in tmp_i2:nstr){
  if(s[i] == "ZCORN"){
    flag1 <- T
    i <- i + 1
  }
  else if(s[i] == "/" & flag1){
    flag1 <- F
    tmp_i3 <- i + 1
    break()
  }
  else if(flag1){
    tmp_row <- as.numeric(strsplit(s[i], split = "\\s+")[[1]])
    
    if(length(tmp_row) >= 1 & !is.na(tmp_row[1])){
      for(i_s in tmp_row){
        z_coords[i_x, i_y, i_z] <- i_s
        i_x <- i_x + 1
        if(i_x > 2*nx){
          i_x <- 1
          i_y <- i_y + 1
          if(i_y > 2*ny){
            i_y <- 1
            i_z <- i_z + 1
            if(i_z > 2*nz) print("OK")
          }
        }
      }
    }
  }
}


##### ACTNUM #####

i_x <- i_y <- i_z <- 1
flag1 <- F

for(i in tmp_i3:nstr){
  if(s[i] == "ACTNUM"){
    flag1 <- T
    i <- i + 1
  }
  else if(s[i] == "/" & flag1){
    flag1 <- F
    tmp_i4 <- i + 1
    break()
  }
  else if(flag1){
    tmp_row <- as.numeric(strsplit(s[i], split = "\\s+")[[1]])
    
    if(length(tmp_row) >= 1 & !is.na(tmp_row[1])){
      for(i_s in tmp_row){
        m_grd[i_x, i_y, i_z] <- i_s
        i_x <- i_x + 1
        if(i_x > nx){
          i_x <- 1
          i_y <- i_y + 1
          if(i_y > ny){
            i_y <- 1
            i_z <- i_z + 1
            if(i_z > nz) print("OK")
          }
        }
      }
    }
  }
}

##### NTG #####

t <- file(paste0(gr_path,gr_name), open = "r")
s <- readLines(t)
close(t)

nstr = length(s)

i_x <- i_y <- i_z <- 1
flag1 <- F

for(i in 1:nstr){
  if(strsplit(s[i], split = "\\s+")[1] == "NTG"){
    flag1 <- T
    i <- i + 1
  }
  else   if(s[i] == "/" & flag1){
    flag1 <- F
    tmp_i4 <- i + 1
    break()
  }
  else if(flag1 & strsplit(s[i], split = "\\s+")[1] != "--"){
    tmp_row <- as.numeric(strsplit(s[i], split = "\\s+")[[1]])
    
    if(length(tmp_row) >= 1 & !is.na(tmp_row[1])){
      for(i_s in tmp_row){
        if(is.numeric(i_s)){
          m_grd[i_x, i_y, i_z] <- i_s
          i_x <- i_x + 1
          if(i_s != 1 & i_s != 0) print(paste(i_s,"WWW"))
          if(i_x > nx){
            i_x <- 1
            i_y <- i_y + 1
            if(i_y > ny){
              i_y <- 1
              i_z <- i_z + 1
              if(i_z > nz) print("OK")
            }
          }
        }
      }
    }
  }
}

##### Find connected volumes #####

options(expressions = 5e5) #Need For large file (standard 5000)

##### Find connected cells at the x_row #####

line_vol <- function(j_x, j_y, j_z, n){
  for(i in j_x:nx){
    if(m_grd[i,j_y,j_z] > 0 &
       m_grd[i,j_y,j_z] == l_grd[i,j_y,j_z]){
      l_grd[i,j_y,j_z] <- n + 1
    }
    else break
  }
  l_grd[j_x,j_y,j_z] <- m_grd[j_x,j_y,j_z]
  for(i in seq(j_x,1,by = -1)){
    if(m_grd[i,j_y,j_z] > 0 &
       m_grd[i,j_y,j_z] == l_grd[i,j_y,j_z]){
      l_grd[i,j_y,j_z] <- n + 1
    }
    else break
  }
}

##### Find connected cells at the xy_plane #####

flat_vol <- function(j_x, j_y, j_z, n){

  line_vol(j_x, j_y, j_z, n)
  
  if(j_y > 1){
    for(i in 1:nx){
      if(l_grd[i,j_y,j_z] == n + 1 &
         m_grd[i,j_y - 1,j_z] > 0 &
         m_grd[i,j_y - 1,j_z] == l_grd[i,j_y - 1,j_z]){
        flat_vol(i, j_y - 1, j_z, n)
      }
    }
  }
  if(j_y < ny){
    for(i in 1:nx){
      if(l_grd[i,j_y,j_z] == n + 1 &
         m_grd[i,j_y + 1,j_z] > 0 &
         m_grd[i,j_y + 1,j_z] == l_grd[i,j_y + 1,j_z]){
        flat_vol(i, j_y + 1, j_z, n)
      }
    }
  }
}

##### Find connected cells at the 3D volume #####

vol_find <- function(j_x, j_y, j_z, n){
  
  flat_vol(j_x, j_y, j_z, n)
  
  if(j_z > 1){
    for(i in j_x:nx){
      for(j in j_y:ny){
        if(l_grd[i,j,j_z] == n + 1 &
           m_grd[i,j,j_z - 1] > 0 &
           m_grd[i,j,j_z - 1] == l_grd[i,j,j_z - 1]){
          vol_find(i, j, j_z - 1, n)
        }
      }
    }
  }
  if(j_z < nz){
    for(i in j_x:nx){
      for(j in j_y:ny){
        if(l_grd[i,j,j_z] == n + 1 &
           m_grd[i,j,j_z + 1] > 0 &
           m_grd[i,j,j_z + 1] == l_grd[i,j,j_z + 1]){
          vol_find(i, j, j_z + 1, n)
        }
      }
    }
  }
}




n_vol <- 0
l_grd <- m_grd

for(i_z in 1:nz){
  for(i_y in 1:ny){
    for(i_x in 1:nx){
      if(m_grd[i_x,i_y,i_z] > 0 &
         m_grd[i_x,i_y,i_z] == l_grd[i_x,i_y,i_z]){
        n_vol <- n_vol + 1
        print(paste(n_vol,i_x,i_y,i_z))
        vol_find(i_x, i_y, i_z, n_vol)
        nnn <- 0
      }
    }
  }
}


##### Prepare STL file #####


add_vert <- function(vert1, vert2, vert3, vert4){
  m_stl[nnn] <<- "  facet normal 0.0 0.0 0.0"
  nnn <<- nnn + 1
  m_stl[nnn] <<- "    outer loop"
  nnn <<- nnn + 1
  
  m_stl[nnn] <<- paste("      vertex", vert1[1], vert1[2], vert1[3])
  nnn <<- nnn + 1
  m_stl[nnn] <<- paste("      vertex", vert2[1], vert2[2], vert2[3])
  nnn <<- nnn + 1
  m_stl[nnn] <<- paste("      vertex", vert3[1], vert3[2], vert3[3])
  nnn <<- nnn + 1
  
  m_stl[nnn] <<- "    endloop"
  nnn <<- nnn + 1
  m_stl[nnn] <<- "  endfacet"
  nnn <<- nnn + 1
  
  m_stl[nnn] <<- "  facet normal 0.0 0.0 0.0"
  nnn <<- nnn + 1
  m_stl[nnn] <<- "    outer loop"
  nnn <<- nnn + 1
  
  m_stl[nnn] <<- paste("      vertex", vert2[1], vert2[2], vert2[3])
  nnn <<- nnn + 1
  m_stl[nnn] <<- paste("      vertex", vert3[1], vert3[2], vert3[3])
  nnn <<- nnn + 1
  m_stl[nnn] <<- paste("      vertex", vert4[1], vert4[2], vert4[3])
  nnn <<- nnn + 1
  
  m_stl[nnn] <<- "    endloop"
  nnn <<- nnn + 1
  m_stl[nnn] <<- "  endfacet"
  nnn <<- nnn + 1
}


add_vert_xyz <- function(j_x, j_y, j_z, z1, z2, z3, z4, side){
  
  if(side == 1){ #low x
    v_x <- c(j_x,
             j_x,
             j_x,
             j_x)
    v_y <- c(j_y,
             j_y + 1,
             j_y,
             j_y + 1)
  }
  else if(side == 2){ #big x
    v_x <- c(j_x + 1,
             j_x + 1,
             j_x + 1,
             j_x + 1)
    v_y <- c(j_y,
             j_y + 1,
             j_y,
             j_y + 1)
  }
  else if(side == 3){ #low y
    v_x <- c(j_x,
             j_x + 1,
             j_x,
             j_x + 1)
    v_y <- c(j_y,
             j_y,
             j_y,
             j_y)
  }
  else if(side == 4){ #big y
    v_x <- c(j_x,
             j_x + 1,
             j_x,
             j_x + 1)
    v_y <- c(j_y + 1,
             j_y + 1,
             j_y + 1,
             j_y + 1)
  }
  else if(side == 5){ #low z
    v_x <- c(j_x,
             j_x,
             j_x + 1,
             j_x + 1)
    v_y <- c(j_y,
             j_y + 1,
             j_y,
             j_y + 1)
  }
  else if(side == 6){ #big z
    v_x <- c(j_x,
             j_x,
             j_x + 1,
             j_x + 1)
    v_y <- c(j_y,
             j_y + 1,
             j_y,
             j_y + 1)
  }
    
  v_z <- c(z1,
           z2,
           z3,
           z4)  
  vert <- list(c(),c(),c(),c())
  
  for(i in 1:4){
    vert[[i]] <- c((pills[v_x[i],v_y[i],1,1] - pills[v_x[i],v_y[i],2,1])*
                   (pills[v_x[i],v_y[i],2,3] - v_z[i])/(pills[v_x[i],v_y[i],2,3] - pills[v_x[i],v_y[i],1,3])+
                   pills[v_x[i],v_y[i],2,1],
                 (pills[v_x[i],v_y[i],1,2] - pills[v_x[i],v_y[i],2,2])*
                   (pills[v_x[i],v_y[i],2,3] - v_z[i])/(pills[v_x[i],v_y[i],2,3] - pills[v_x[i],v_y[i],1,3])+
                   pills[v_x[i],v_y[i],2,2],
                 max_depth - v_z[i])
  }
  
  add_vert(vert[[1]],vert[[2]],vert[[3]],vert[[4]])
}

# Coords norm

x_min <- min(pills[,,,1])
x_max <- max(pills[,,,1])
y_min <- min(pills[,,,2])
y_max <- max(pills[,,,2])
z_min <- min(pills[,,,3])
z_max <- max(pills[,,,3])

# Invert Y
pills[,,,2] <- y_max - pills[,,,2]

L_max <- max(x_max-x_min,y_max-y_min)

x_max <- x_min + L_max
y_max <- y_min + L_max

pills[,,,1] <- (pills[,,,1] - x_min) * 1000 / x_max
pills[,,,2] <- (pills[,,,2] - y_min) * 1000 / y_max
pills[,,,3] <- (pills[,,,3] - z_min) * 1000 / z_max


library(data.table)


for(i in 2:max(l_grd)){
  
  m_stl <- list()
  nnn <- 1
  
  m_stl[nnn] <- paste0("solid layer_",i-1)
  nnn <- nnn + 1
  for(i_z in 1:nz){
    print(paste0("Object: ",i-1,"/",max(l_grd)-1," Z: ", round(i_z/nz*100, digits = 1),"%"))
    for(i_y in 1:ny){
      for(i_x in 1:nx){
        if(l_grd[i_x,i_y,i_z] == i){
          if(i_x == 1){
            z1 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z - 1]
            z2 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z - 1]
            z3 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z]
            z4 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,1)
          }
          else if(l_grd[i_x - 1,i_y,i_z] == 0){
            z1 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z - 1]
            z2 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z - 1]
            z3 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z]
            z4 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,1)
          }
          if(i_x == nx){
            z1 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z - 1]
            z2 <- z_coords[2*i_x, 2*i_y, 2*i_z - 1]
            z3 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z]
            z4 <- z_coords[2*i_x, 2*i_y, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,2)
          }
          else if(l_grd[i_x + 1,i_y,i_z] == 0){
            z1 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z - 1]
            z2 <- z_coords[2*i_x, 2*i_y, 2*i_z - 1]
            z3 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z]
            z4 <- z_coords[2*i_x, 2*i_y, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,2)
          }
          
          if(i_y == 1){
            z1 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z - 1]
            z2 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z - 1]
            z3 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z]
            z4 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,3)
          }
          else if(l_grd[i_x,i_y - 1,i_z] == 0){
            z1 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z - 1]
            z2 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z - 1]
            z3 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z]
            z4 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,3)
          }
          if(i_y == ny){
            z1 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z - 1]
            z2 <- z_coords[2*i_x, 2*i_y, 2*i_z - 1]
            z3 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z]
            z4 <- z_coords[2*i_x, 2*i_y, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,4)
          }
          else if(l_grd[i_x,i_y + 1,i_z] == 0){
            z1 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z - 1]
            z2 <- z_coords[2*i_x, 2*i_y, 2*i_z - 1]
            z3 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z]
            z4 <- z_coords[2*i_x, 2*i_y, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,4)
          }
          
          if(i_z == 1){
            z1 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z - 1]
            z2 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z - 1]
            z3 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z - 1]
            z4 <- z_coords[2*i_x, 2*i_y, 2*i_z - 1]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,5)
          }
          else if(l_grd[i_x,i_y,i_z - 1] == 0){
            z1 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z - 1]
            z2 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z - 1]
            z3 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z - 1]
            z4 <- z_coords[2*i_x, 2*i_y, 2*i_z - 1]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,5)
          }
          if(i_z == nz){
            z1 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z]
            z2 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z]
            z3 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z]
            z4 <- z_coords[2*i_x, 2*i_y, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,6)
          }
          else if(l_grd[i_x,i_y,i_z + 1] == 0){
            z1 <- z_coords[2*i_x - 1, 2*i_y - 1, 2*i_z]
            z2 <- z_coords[2*i_x - 1, 2*i_y, 2*i_z]
            z3 <- z_coords[2*i_x, 2*i_y - 1, 2*i_z]
            z4 <- z_coords[2*i_x, 2*i_y, 2*i_z]
            add_vert_xyz(i_x,i_y,i_z,z1,z2,z3,z4,6)
          }
        }
      }
    }
  }
  m_stl[nnn] <- "endsolid"
  nnn <- nnn + 1
  
  fwrite(m_stl, file = paste0(gr_path,"M_",i-1,".stl"), sep = "\n")
}