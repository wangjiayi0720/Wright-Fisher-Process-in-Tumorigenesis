# parameters
n0<- 10^6
n1<- 10^9
a<- 1.72*10^(-3)
u<- 10^(-7)
d<- 100
s<- 0.01
C<- (log(s/u/d))^2/s/log(n0*n1)
# simulation times
simu<- 50

# time recorded for each mutation first presents
time<- matrix(0,simu,20)

# transition probability p
p<- rep(0,101)
# weight w
w<- NULL
for (k in 1:101)
{
  w[k]=(1+s)^(k-1)
}


for(m in 1:simu)
{
  # number of cells in each mutation group
  cell<- matrix(0,1,101)
  cell[1,1]<- n0
  temp<- cell
  # time t
  t<- 0
  for(tao in 1:20)
  {
    if(tao <=10)
    {
      while(temp[(tao+1)]==0)
      {
        t<- t+1
        # calculate each time transition probability
        for (k in 0:100)
        {
          x<- 0
          for(j in 0:k)
          {
            x<- x + dbinom(k-j,d-j,u)*w[j+1]*temp[j+1]/sum(w*temp)
          }
          p[k+1]<- x
        }
        temp<- t(rmultinom(1,n0*(1+a)^t,p))
      }
      time[m,tao]<- t
    }
    else
    {
      while(temp[(tao+1)]==0)
      {
        t<- t+1
        # calculate each time transition probability
        for (k in 0:100)
        {
          x<- 0
          for(j in 0:k)
          {
            x<- x + dbinom(k-j,d-j,u)*w[j+1]*temp[j+1]/sum(w*temp)
          }
          p[k+1]<- x
        }   
        n_total <- n0*(1+a)^t
        n_100 <- n_total-n_total%%500
        n_chunk <- n_100 / 500
        temp<- rep(0,101)
        for(i in 1:500)
        {
          temp <- temp + t(rmultinom(1,n_chunk,p))
        }
        temp<- temp+t(rmultinom(1,n_total%%500,p))
      }
      time[m,tao]<- t
    }                      
  }
}

