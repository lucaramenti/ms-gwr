#The required libraries are loaded----
library(plot3D)
library(GWmodel)
library(psych)
library(ggplot2)
library(reshape2)

#The required functions are loaded----
source("functions.R")
#We build two grids in a range from -5 to 5 with step equal to 0.5, in order to have station- and event-coordinates----
inf = -5
sup = 5
step = 0.5
x = s1 = s2 = e1 = e2 = seq(inf, sup, by = step)
adj = (length(x)+1)/2

grid_e_sim = matrix(0,(length(e1)*length(e2)),2)
for (i in 1:length(e1)){
  for (j in 1:length(e2)){
    grid_e_sim[j+(i-1)*length(s2),] = c(e1[i],e2[j])
  }
}

grid_s_sim = matrix(0,(length(s1)*length(s2)),2)
for (i in 1:length(s1)){
  for (j in 1:length(s2)){
    grid_s_sim[j+(i-1)*length(s2),] = c(s1[i],s2[j])
  }
}

#The regression coefficients are set and, in the case of spatial non-stationarity, they are also plotted----
#beta_c intercept
beta_c_intercept_true = 8

#beta_c variable
beta_c_variable_true = 4

x11()
par(mfrow=c(1,2), mai = c(0.4, 0.4, 0.4, 0.4))

#beta_e variable
beta_e_variable_true = matrix(0,length(e1), length(e2))
beta_e_variable_true_col = matrix(0,length(e1)*length(e2),3)
for (i in 1:length(e1)){
  for (j in 1:length(e2)){
    #beta_e_variable_true[i,j] = -0.01*(s1[i]+0.5)^3 - 0.01*(s2[j]-0.5)^3 + 0.007*(s2[j]-2.5)^3+3.5
    beta_e_variable_true[i,j] = +0.05*(e1[i]+1)^2 + 0.05*(e2[j])^2 + 1.5 + 0.1*(e1[i]-2)
    beta_e_variable_true_col[j+(i-1)*length(e2),] = c(beta_e_variable_true[i,j],e1[i],e2[j])
  }
}

persp3D(e1,e2,beta_e_variable_true, colvar = beta_e_variable_true,theta=20,
        phi=20,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        #contour = list(col = "grey", side = c("zmin", "z")),
        main=expression(paste(beta["1E"]))
)

#beta_s variable
beta_s_variable_true = matrix(0,length(s1), length(s2))
beta_s_variable_true_col = matrix(0,length(s1)*length(s2),3)
for (i in 1:length(s1)){
  for (j in 1:length(s2)){
    #beta_s_variable_true[i,j] = -0.05*(s1[i]-1)^2 - 0.05*(s2[j])^2 + 4.5 - 0.1*(s1[i]+3)
    beta_s_variable_true[i,j] = -0.01*(s1[i]+0.5)^3 - 0.01*(s2[j]-0.5)^3 + 0.007*(s2[j]-2.5)^3+3.5
    beta_s_variable_true_col[j+(i-1)*length(s2),] = c(beta_s_variable_true[i,j],s1[i],s2[j])
  }
}

persp3D(s1,s2,beta_s_variable_true, colvar = beta_s_variable_true,theta=20,
        phi=20,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        #contour = list(col = "grey", side = c("zmin", "z")),
        main=expression(paste(beta["1S"]))
)

#Now we have to create the dataset, along some normal distributions, considering a stationary intercept----
#and three covariates, a stationary one, a site-dependent one and an event-dependent one. The response variables 
#are also computed, after introducing the error term.
#create X
set.seed(6494)
n = length(s1)^4
Xc_var = rnorm(n, 3, 5)
Xe_var = rnorm(n, -4, 5)
Xs_var = rnorm(n, 4, 5)

Xc = cbind(rep(1,n), Xc_var)
Xe = Xe_var
Xs = Xs_var

#create variance on four dimensional grid
epsilon = rnorm(n,0,2) 

n = length(s1)^4 #create "four-dimensional grid"
f = 4
l = rep(list(x), f)
coords = expand.grid(l)
X_coords = cbind(Xc_var, Xe_var, Xs_var, coords, epsilon)
names(X_coords)[4:7]=c("E1", "E2", "S1", "S2")

#response y
beta_c_true = cbind(beta_c_intercept_true, beta_c_variable_true)
beta_s_true = beta_s_variable_true_col
beta_e_true = beta_e_variable_true_col
yc = rep(0,n)
for (i in 1:n){
  coords_e = X_coords[i,4:5]
  coords_s = X_coords[i,6:7]
  yc[i] = Xc[i,]%*%t(beta_c_true) +
    Xe_var[i]*beta_e_variable_true[coords_e[[1]]/step+adj,coords_e[[2]]/step+adj] +
    Xs_var[i]*beta_s_variable_true[coords_s[[1]]/step+adj,coords_s[[2]]/step+adj] +
    epsilon[i]
}


#The next step is to create our simulation subset, made of 80 elements picked casually from our generated dataset----
set.seed(2803)
n_sample = 80
lines = sample(1:n, n_sample)
ordered_lines = sort(lines)
X_sim = X_coords[ordered_lines,]
y_sim_c = yc[ordered_lines]

grid_sim = matrix(0,(length(s1)*length(s2)),2)
for (i in 1:length(s1)){
  for (j in 1:length(s2)){
    grid_sim[j+(i-1)*length(s2),] = c(s1[i],s2[j])
  }
}

coords_e_sim = cbind(X_sim$E1, X_sim$E2)
coords_s_sim = cbind(X_sim$S1, X_sim$S2)

#We compute the best bandwidth for the complete model----
#since we assume it to be the correct one; morover we see if it is better 
#to use ESC or SEC. Notice that in this example we force bw_e = bw_s
bw_sec = bw_cv(0.5,3.2,0.15,80,SEC_only_calibration, "sec", X_sim$Xc_var, X_sim$Xe_var, X_sim$Xs_var, y_sim_c,
         "c", coords_e_sim, coords_s_sim)
bw_esc = bw_cv(0.5,3.2,0.15,80,ESC_only_calibration, "esc", X_sim$Xc_var, X_sim$Xe_var, X_sim$Xs_var, y_sim_c,
         "c", coords_e_sim, coords_s_sim)
cvsum_sec = bw_sec[[2]][,2]
cvsum_esc = bw_esc[[2]][,2]
bws = bw_esc[[2]][,1]
dataset = as.data.frame(cbind(bws, cvsum_sec, cvsum_esc))
melted = melt(dataset, id.vars = "bws")

x11()
ggplot(melted, aes(x=bws, y=value, col=variable)) + 
  geom_line(size = 1.7)+
  labs(y = "CVSS", x = "Bandwidth")+
  scale_color_manual(values=c("royalblue3","skyblue"),labels = c("SEC", "ESC"))+
  theme(axis.text=element_text(size=15, colour = "black"),
        axis.title=element_text(size=25),
        legend.text=element_text(size=25, colour = "black"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA)
  )+
  scale_x_continuous(breaks = seq(from = 0.65, to = 3.05, by = 0.3)) +
  geom_vline(xintercept = 1.25, size = 1.2, colour="red3")+
  coord_cartesian(xlim = c(0.6, 3.1),ylim = c(800, 3400))+
  geom_point(size = 4)
bwe = bws = 1.25 #computed with code before, results are saved
#Both ESC and SEC lead to bw_e = bw_s = 1.25, but ESC performs better as far as CVSS is concerned.
#Now we proceed by carrying out the permutation test, to verify if we are considering the right model.
#To carry out all these tests, we use the bandwidth which are found for the complete model in the previous point.
#At first we verify if it is reasonable to introduce spatial non-stationarity----
#OLS vs completely varying
ols = lm(y_sim_c ~ X_sim$Xc_var + X_sim$Xe_var + X_sim$Xs_var)
bwe = bws = 1.25
esc_only_intercept = ESC_only_constant_intercept_calibration(cbind(X_sim$Xc_var, X_sim$Xe_var), X_sim$Xs_var, y_sim_c, bwe, bws, coords_e_sim, coords_s_sim)
#compute R(H0)
X = cbind(rep(1,n_sample), X_sim$Xc_var, X_sim$Xe_var, X_sim$Xs_var)
I = diag(1,n_sample)
Hols = X%*%(solve(t(X)%*%X))%*%t(X)
RH0 = t(I-Hols)%*%(I-Hols)
eps = (I-Hols)%*%y_sim_c
#compute R(H1)
B = esc_only_intercept$B
Xcc = rep(1,dim(X_sim)[1])
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
set.seed(6494)
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = Hols %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
}
p = sum(t_stat>as.numeric(T0))/n_perm
p

#Then we check if it is correct to consider Xc as constant and Xe and Xs as spatially varying----
#check if Xc is constant
#calibrate models under H0 and H1
esc = ESC_only_calibration(X_sim$Xc_var, X_sim$Xe_var, X_sim$Xs_var, y_sim_c, "c", bwe, bws, coords_e_sim, coords_s_sim)
esc_only_intercept = ESC_only_constant_intercept_calibration(cbind(X_sim$Xc_var, X_sim$Xe_var), X_sim$Xs_var, y_sim_c, bwe, bws, coords_e_sim, coords_s_sim)
#compute R(H0)
I= diag(1,n_sample)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
B = esc_only_intercept$B
Xcc = rep(1,dim(X_sim)[1])
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
set.seed(6494)
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
}
p_xc = sum(t_stat>as.numeric(T0))/n_perm
p_xc

#check if Xe is constant
#calibrate models under H0 and H1
esc = ESC_only_calibration(X_sim$Xe_var, X_sim$Xc_var, X_sim$Xs_var, y_sim_c, "c", bwe, bws, coords_e_sim, coords_s_sim)
esc_only_intercept = ESC_only_constant_intercept_calibration(cbind(X_sim$Xe_var, X_sim$Xc_var), X_sim$Xs_var, y_sim_c, bwe, bws, coords_e_sim, coords_s_sim)
#compute R(H0)
I= diag(1,n_sample)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xe_var)
H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
B = esc_only_intercept$B
Xcc = rep(1,dim(X_sim)[1])
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
}
p_xe = sum(t_stat>as.numeric(T0))/n_perm
p_xe

#check if Xs is constant
#calibrate models under H0 and H1
esc = mixed_SC_calibration(X_sim$Xs_var, cbind(X_sim$Xe_var, X_sim$Xc_var), y_sim_c, 1.25, "c", coords_e_sim)
esc_only_intercept = ESC_only_constant_intercept_calibration(cbind(X_sim$Xe_var, X_sim$Xc_var), X_sim$Xs_var, y_sim_c, bwe, bws, coords_e_sim, coords_s_sim)
#compute R(H0)
I= diag(1,n_sample)
Hs = esc$Hs
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xs_var)
# H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
H0 = Hs + (I-Hs)%*%Xcc%*%solve(t(Xcc)%*%t(I-Hs)%*%(I-Hs)%*%Xcc)%*%t(Xcc)%*%t(I-Hs)%*%(I-Hs)
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
B = esc_only_intercept$B
Xcc = rep(1,dim(X_sim)[1])
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
}
p_xs = sum(t_stat>as.numeric(T0))/n_perm
p_xs

#The next step is to check if the constant terms can be considered as null----
#check if the intercept is null
#calibrate models under H0 and H1
esc_null = mixed_SC_no_intercept_calibration(X_sim$Xc_var, X_sim$Xs_var, y_sim_c, 1.25, coords_s_sim)
esc = ESC_only_calibration(X_sim$Xc_var, X_sim$Xe_var, X_sim$Xs_var, y_sim_c, "c", bwe, bws, coords_e_sim, coords_s_sim)
#compute R(H0)
I= diag(1,n_sample)
Hs = esc_null$Hs
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H0 = Hs + (I-Hs)%*%Xcc%*%solve(t(Xcc)%*%t(I-Hs)%*%(I-Hs)%*%Xcc)%*%t(Xcc)%*%t(I-Hs)%*%(I-Hs)
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
I= diag(1,n_sample)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
}
p_intercept_null = sum(t_stat>as.numeric(T0))/n_perm
p_intercept_null

#check if Xc is null
#calibrate models under H0 and H1
esc_null = mixed_SC_no_intercept_calibration(rep(1,length(X_sim$Xe_var)), X_sim$Xe_var, y_sim_c, 1.25, coords_e_sim)
esc = ESC_only_calibration(X_sim$Xc_var, X_sim$Xe_var, X_sim$Xs_var, y_sim_c, "c", bwe, bws, coords_e_sim, coords_s_sim)

#compute R(H0)
I= diag(1,n_sample)
Hs = esc_null$Hs
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H0 = Hs + (I-Hs)%*%Xcc%*%solve(t(Xcc)%*%t(I-Hs)%*%(I-Hs)%*%Xcc)%*%t(Xcc)%*%t(I-Hs)%*%(I-Hs)
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
I= diag(1,n_sample)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
}
p_xc_null = sum(t_stat>as.numeric(T0))/n_perm
p_xc_null

#We can now execute the whole calibration, using the previously obtained results to set parameters and methods----
bwe = bws = 1.25
esc = ESC_general(X_sim$Xc_var, X_sim$Xe_var, X_sim$Xs_var, y_sim_c, "c", bwe, bws, coords_e_sim, coords_s_sim, as.matrix(grid_sim))

betas_c = esc$beta_c
betas_c

#The spatially varying coefficients can be graphically compared with the true ones----
beta_s = esc$beta_s
beta_e = esc$beta_e

beta_s_grid = matrix(0, length(s1),length(s2))
k=1
for (i in 1:length(s1)){
  for (j in 1:length(s2)){
    beta_s_grid[i,j] = beta_s[1,k]
    k = k+1
  }
}

beta_e_grid = matrix(0, length(e1),length(e2))
k=1
for (i in 1:length(e1)){
  for (j in 1:length(e2)){
    beta_e_grid[i,j] = beta_e[1,k]
    k = k+1
  }
}

x11()
par(mfrow=c(2,2),  mai = c(0.4, 0.4, 0.4, 0.4))
persp3D(e1,e2,beta_e_variable_true, colvar = beta_e_variable_true,theta=30,
        phi=40,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        main=expression(paste(beta["1E"]))
)

persp3D(s1,s2,beta_s_variable_true, colvar = beta_s_variable_true,theta=30,
        phi=40,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        main=expression(paste(beta["1S"]))
)

persp3D(e1,e2,beta_e_grid, colvar = beta_e_grid,theta=30,
        phi=40,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        main=expression(paste(beta["1E, EMGWR"]))
)

persp3D(s1,s2,beta_s_grid, colvar = beta_s_grid,theta=30,
        phi=40,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        main=expression(paste(beta["1S, EMGWR"]))
)

#Finally we can compute R^2_{adj}----
I= diag(1,n_sample)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H_esc = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
res = (I-H_esc)%*%y_sim_c
delta1 = n_sample-2*tr(H_esc)+tr(t(H_esc)%*%H_esc)
rss = sum(res^2)
sigma2hat = rss/delta1
tss = sum((H_esc%*%y_sim_c-mean(y_sim_c))^2)
R2 = 1-rss/tss
R2adj = 1-(1-R2)*(n_sample-1)/delta1
R2adj