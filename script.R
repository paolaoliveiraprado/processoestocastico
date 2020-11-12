# PAOLA DE OLIVEIRA PRADO

# x0 = passo inicial
# n  = até onde o processo vai
# k  = alcance
# p  = matriz de transição
# r  = número de replicações

# Função para simulação de uma cadeira de alcance 0, 1, 2 e 3.

tarefa=function(x0,n,k,p,r){
  aux=c(x0)
  xn=c(0)
  rep=matrix(c(0), nrow=r, ncol=(n), byrow=TRUE)
  for(j in 1:r){ # vai preencher a linha de cada repetição
    u=runif(n)
    for(i in 1:n){
      if(k==0){ #PARA k==0
        if(u[i]<p[1,1]){
          xn[i]=0
        } else { 
          xn[i]=1
        }
        aux=xn[i]
      } 
      else if (k==1){   #PARA k==1
        if(aux==0){
          if(u[i]<p[1,1]){
            xn[i]=0
          } else {
            xn[i]=1
          }
        }
        else{
          if(u[i]<p[2,1]){
            xn[i]=0
          } else {
            xn[i]=1
          }
        }
        aux=xn[i]
      } 
      else if (k==2){   #PARA k==2
        if(aux[1]==0 & aux[2]==0){
          if(u[i]<p[1,1]){
            xn[i]=0
          } else {
            xn[i]=1
          }
        }
        else if (aux[1]==0 & aux[2]==1){
          if(u[i]<p[2,1]){
            xn[i]=0
          } else {
            xn[i]=1
          }
        }
        else if(aux[1]==1 & aux[2]==0){
          if(u[i]<p[3,1]){
            xn[i]=0
          } else {
            xn[i]=1
          }        
        }
        else {
          if(u[i]<p[4,1]){
            xn[i]=0
          } else {
            xn[i]=1
          }        
        }
        aux=c(aux[2],xn[i])
      }
      else { #PARA K==3
        if(aux[1]==0 & aux[2]==0 & aux[3]==0){
          if(u[i]<p[1,1]){
            xn[i]=0 
          }else {
            xn[i]=1
          }
        }
        else if(aux[1]==0 & aux[2]==0 & aux[3]==1){
          if(u[i]<p[2,1]){
            xn[i]=0
          } else {
            xn[i]=1
          } 
        }
        else if(aux[1]==0 & aux[2]==1 & aux[3]==0){
          if(u[i]<p[3,1]){
            xn[i]=0
          } else{
            xn[i]=1
          }
        }
        else if(aux[1]==0 & aux[2]==1 & aux[3]==1){
          if(u[i]<p[4,1]){
            xn[i]=0
          } else{
            xn[i]=1
          }
        }
        else if(aux[1]==1 & aux[2]==0 & aux[3]==0){
          if(u[i]<p[5,1]){
            xn[i]=0
          } else{
            xn[i]=1
          }
        }
        else if(aux[1]==1 & aux[2]==0 & aux[3]==1){
          if(u[i]<p[6,1]){
            xn[i]=0
          } else {
            xn[i]=1
          }
        }
        else if(aux[1]==1 & aux[2]==1 & aux[3]==0){
          if(u[i]<p[7,1]){
            xn[i]=0
          } else {
            xn[i]=1
          }
        }
        else{
          if(u[i]<p[8,1]){
            xn[i]=0
          } else{
            xn[i]=1
          }
        }
        aux=c(aux[2],aux[3],xn[i])
        }
      }
    rep[j,]=matrix(c(xn),nrow=1, ncol=(n), byrow=TRUE) #tem que estar dps de terminar o for i 
  }
  return(rep) #lista de replicações cada linha(n) é uma replicação
}

r=100      #replicações- fixo!
n=100      #tamanho da amostra, 100, 1000 ou 10000

#####  RETIRAR O # DA MATRIZ DE TRANSIÇÃO (p) E DA SIMULAÇÃO DE ACORDO COM O k QUE DESEJA SIMULAR ! ##### 

##### Simulando cadeia para k=0 #####

#p=matrix(c(0.3,0.7), ncol=2, byrow=T);p #MATRIZ DE TRANSIÇÃO para k=0
#k=tarefa(x0=NA,n,k=0,p,r)        #Simulação para k=0


##### Simulando cadeia para k=1 #####

#p=matrix(c(0.4,0.6,0.2,0.8), ncol=2, byrow=T) #MATRIZ DE TRANSIÇÃO PARA k=1
#k=tarefa(x0=1,n,k=1,p,r)         #Simulação para k=1


##### Simulando cadeia para k=2 #####

#p=matrix(c(0.5,0.5,0.2,0.8,0.3,0.7,0.4,0.6), ncol=2 , byrow = T);p #Matriz de transição para k=2
#k=tarefa(x0=c(1,1),n,k=2,p,r)    #Simulação para k2

##### Simulando cadeia para k=3 #####


k=tarefa(x0=c(0,0,0),n,k=3,p,r)  #Simulação para k=3
p=matrix(c(0.5,0.5,0.4,0.6,0.4,0.6,0.3,0.7,0.5,0.5,0.3,0.7,0.3,0.7,0.3,0.7), ncol=2,byrow=T);p #Matriz de transição para k=3


k #Matriz com a evolução temporal de todas as replicações

##### 3) Implementando o BIC - para K=0 ; n=100 ; n=1000 ; n=10000 #####

# Fórmula DO BIC , para K=0 : BIC_0 = prob0^n0 * prob1^n1 - 1/2 (2^0 * (abs(2)-1) * log(n))

N0=c(0)
prob0=c(0)
N1=c(0)
prob1=c(0)
ver0=c(0)
for(i in 1:r){
  N0[i]=sum(k[i,] == 0)
  prob0[i]=N0[i]/n
  N1[i]=sum(k[i,] == 1)
  prob1[i]=N1[i]/n
  ver0[i]=log2((prob0[i]^N0[i]) * (prob1[i]^N1[i]))
}
N0
N1
prob0 #probabilidade em cada replicação de sair o valor 0
prob1

##### BIC k=0 ####

bic0=c(0)
for(i in 1:r){
  bic0[i]= ver0[i] - 1/2 * (2^0 * (abs(2)-1) * log(n,2))
}
bic0

##### 3) Implementando o BIC - para K=1 ; n=100 ; n=1000 ; n=10000 #####
# Fórmula DO BIC , para K=1 : BIC_1 = prob(0|0)^n00 * prob(1|0)^n01 * prob(0|1)^ * prob1 - 1/2 (2^1 * (abs(2)-1) * log(n))


#Contagem de N00 para cada replicação
N00 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-1)){  #1:R-1
    if(k[i,j]==0 & k[i,j+1] == 0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N00[i]=cont
  cont=0
}
N00

#Contagem de N01 para cada replicação
N01 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-1)){  #1:R-1
    if(k[i,j]==0 & k[i,j+1] == 1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N01[i]=cont
  cont=0
}
N01

#Contagem de N10 para cada replicação
N10 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-1)){  #1:R-1
    if(k[i,j]==1 & k[i,j+1] == 0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N10[i]=cont
  cont=0
}
N10

#Contagem de N11 para cada replicação
N11 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-1)){  #1:R-1
    if(k[i,j]==1 & k[i,j+1] == 1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N11[i]=cont
  cont=0
}
N11


## Probabidade de sair 00

p00=c(0)
for(i in 1:r){
  p00[i]= N00[i] / (N00[i] + N01[i])
}
p00


## Probabidade de sair 01

p01=c(0)
for(i in 1:r){
  p01[i]= N01[i] / (N01[i] + N00[i])
}
p01


## Probabidade de sair 10

p10=c(0)
for(i in 1:r){
  p10[i]= N10[i] / (N10[i] + N11[i])
}
p10

## Probabidade de sair 11

p11=c(0)
for(i in 1:r){
  p11[i]= N11[i] / (N11[i] + N10[i])
}
p11

##### Verossimilhança k=1 ####

ver1=c(0)
for(i in 1:r){
  ver1[i]= log2(p00[i]^N00[i] * p01[i]^N01[i] * p11[i]^N11[i] * p10[i]^N10[i])
}
ver1

##### BIC k=1 ####

bic1=c(0)
for(i in 1:r){
  bic1[i]= ver1[i] - 1/2 * (2^1 * (abs(2)-1) * log(n,2)) 
}
bic1



##### 3) Implementando o BIC - para K=2 ; n=100 ; n=1000 ; n=10000 #####
# Fórmula DO BIC , para K=2 : BIC_2 = ..... - 1/2 (2^2 * (abs(2)-1) * log(n))

#Contagem de N000 para cada replicação
N000 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-2)){  #1:R-2
    if(k[i,j]==0 & k[i,j+1] == 0 & k[i,j+2]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N000[i]=cont
  cont=0
}
N000

#Contagem de N001 para cada replicação
N001 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-2)){  #1:R-2
    if(k[i,j]==0 & k[i,j+1] == 0 & k[i,j+2]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N001[i]=cont
  cont=0
}
N001

#Contagem de N100 para cada replicação
N100 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-2)){  #1:R-2
    if(k[i,j]==1 & k[i,j+1] == 0 & k[i,j+2]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N100[i]=cont
  cont=0
}
N100

#Contagem de N101 para cada replicação
N101 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-2)){  #1:R-2
    if(k[i,j]==1 & k[i,j+1] == 0 & k[i,j+2]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N101[i]=cont
  cont=0
}
N101

#Contagem de N010 para cada replicação
N010 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-2)){  #1:R-2
    if(k[i,j]==0 & k[i,j+1] == 1 & k[i,j+2]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N010[i]=cont
  cont=0
}
N010

#Contagem de N011 para cada replicação
N011 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-2)){  #1:R-2
    if(k[i,j]==0 & k[i,j+1] == 1 & k[i,j+2]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N011[i]=cont
  cont=0
}
N011

#Contagem de N110 para cada replicação
N110 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-2)){  #1:R-2
    if(k[i,j]==1 & k[i,j+1] == 1 & k[i,j+2]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N110[i]=cont
  cont=0
}
N110

#Contagem de N111 para cada replicação
N111 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-2)){  #1:R-2
    if(k[i,j]==1 & k[i,j+1] == 1 & k[i,j+2]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N111[i]=cont
  cont=0
}
N111

## Probabidade de sair 000 e 001

p000=c(0)
p001=c(0)
for(i in 1:r){
  p000[i]= N000[i] / (N000[i] + N001[i])
  p001[i]=abs(p000[i]-1)
}
p000
p001

## Probabidade de sair 100 e 101

p100=c(0)
p101=c(0)
for(i in 1:r){
  p100[i]= N100[i] / (N100[i] + N101[i])
  p101[i]=abs(p100[i]-1)
}

## Probabidade de sair 010 e 011

p010=c(0)
p011=c(0)
for(i in 1:r){
  p010[i]= N010[i] / (N010[i] + N011[i])
  p011[i]=abs(p010[i]-1)
}

## Probabidade de sair 010 e 011

p010=c(0)
p011=c(0)
for(i in 1:r){
  p010[i]= N010[i] / (N010[i] + N011[i])
  p011[i]=abs(p010[i]-1)
}

## Probabidade de sair 110 e 111

p110=c(0)
p111=c(0)
for(i in 1:r){
  p110[i]= N110[i] / (N110[i] + N111[i])
  p111[i]=abs(p110[i]-1)
}

##### Verossimilhança k=2 ####

ver2=c(0)
for(i in 1:r){
  ver2[i]= log2(p000[i]^N000[i] * p001[i]^N001[i] * p100[i]^N100[i] * p101[i]^N101[i] * p010[i]^N010[i] * p011[i]^N011[i] * p110[i]^N110[i] * p111[i]^N111[i]) 
}
#vero = 8 multiplicações 
ver2

##### BIC k=2 ####

bic2=c(0)
for(i in 1:r){
  bic2[i]= ver2[i] - 1/2 * (2^2 * (abs(2)-1) * log(n,2))
}
bic2

##### 3) Implementando o BIC - para K=3 ; n=100 ; n=1000 ; n=10000 #####
# Fórmula DO BIC , para K=3 : BIC_3 = ..... - 1/2 (2^3 * (abs(2)-1) * log(n)) #16 contagens ! 

#Contagem de N0000 para cada replicação !!!!

N0000 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==0 & k[i,j+1] == 0 & k[i,j+2]==0 & k[i,j+3]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N0000[i]=cont
  cont=0
}
N0000

#Contagem de N0001 para cada replicação

N0001 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==0 & k[i,j+1] == 0 & k[i,j+2]==0 & k[i,j+3]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N0001[i]=cont
  cont=0
}
N0001

#Contagem de N0010 para cada replicação

N0010 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==0 & k[i,j+1] == 0 & k[i,j+2]==1 & k[i,j+3]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N0010[i]=cont
  cont=0
}
N0010

#Contagem de N0011 para cada replicação

N0011 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==0 & k[i,j+1] == 0 & k[i,j+2]==1 & k[i,j+3]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N0011[i]=cont
  cont=0
}
N0011

#Contagem de N0100 para cada replicação

N0100 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==0 & k[i,j+1]==1 & k[i,j+2]==0 & k[i,j+3]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N0100[i]=cont
  cont=0
}
N0100

#Contagem de N0101 para cada replicação

N0101 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==0 & k[i,j+1] ==1 & k[i,j+2]==0 & k[i,j+3]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N0101[i]=cont
  cont=0
}
N0101

#Contagem de N0110 para cada replicação

N0110 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==0 & k[i,j+1] ==1 & k[i,j+2]==1 & k[i,j+3]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N0110[i]=cont
  cont=0
}
N0110

#Contagem de N0111 para cada replicação

N0111 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==0 & k[i,j+1]==1 & k[i,j+2]==1 & k[i,j+3]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N0111[i]=cont
  cont=0
}
N0111

#Contagem de N1000 para cada replicação

N1000 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==1 & k[i,j+1]==0 & k[i,j+2]==0 & k[i,j+3]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N1000[i]=cont
  cont=0
}
N1000

#Contagem de N1001 para cada replicação

N1001 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==1 & k[i,j+1]==0 & k[i,j+2]==0 & k[i,j+3]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N1001[i]=cont
  cont=0
}
N1001

#Contagem de N1010 para cada replicação

N1010 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==1 & k[i,j+1]==0 & k[i,j+2]==1 & k[i,j+3]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N1010[i]=cont
  cont=0
}
N1010

#Contagem de N1011 para cada replicação

N1011 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==1 & k[i,j+1]==0 & k[i,j+2]==1 & k[i,j+3]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N1011[i]=cont
  cont=0
}
N1011

#Contagem de N1100 para cada replicação

N1100 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==1 & k[i,j+1]==1 & k[i,j+2]==0 & k[i,j+3]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N1100[i]=cont
  cont=0
}
N1100

#Contagem de N1101 para cada replicação

N1101 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==1 & k[i,j+1]==1 & k[i,j+2]==0 & k[i,j+3]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N1101[i]=cont
  cont=0
}
N1101

#Contagem de N1110 para cada replicação

N1110 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==1 & k[i,j+1]==1 & k[i,j+2]==1 & k[i,j+3]==0){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N1110[i]=cont
  cont=0
}
N1110

#Contagem de N1111 para cada replicação

N1111 = c(0)
cont=0
valor=0
for(i in 1:r){ #1:R
  for(j in 1:(r-3)){  #1:R-3
    if(k[i,j]==1 & k[i,j+1]==1 & k[i,j+2]==1 & k[i,j+3]==1){
      valor[i]=cont +1  
      cont=valor[i]
    }
  }
  N1111[i]=cont
  cont=0
}
N1111

## Probabidade de sair 0000 e 0001

p0000=c(0)
p0001=c(0)
for(i in 1:r){
  p0000[i]= N0000[i] / (N0000[i] + N0001[i])
  p0001[i]=abs(p0000[i]-1)
}
p0000
p0001

## Probabidade de sair 0010 e 0011

p0010=c(0)
p0011=c(0)
for(i in 1:r){
  p0010[i]= N0010[i] / (N0010[i] + N0011[i])
  p0011[i]=abs(p0010[i]-1)
}
p0010
p0011

## Probabidade de sair 0100 e 0101

p0100=c(0)
p0101=c(0)
for(i in 1:r){
  p0100[i]= N0100[i] / (N0100[i] + N0101[i])
  p0101[i]=abs(p0100[i]-1)
}
p0010
p0011

## Probabidade de sair 0110 e 0111

p0110=c(0)
p0111=c(0)
for(i in 1:r){
  p0110[i]= N0110[i] / (N0110[i] + N0111[i])
  p0111[i]=abs(p0110[i]-1)
}
p0110
p0111

## Probabidade de sair 1000 e 1001

p1000=c(0)
p1001=c(0)
for(i in 1:r){
  p1000[i]= N1000[i] / (N1000[i] + N1001[i])
  p1001[i]=abs(p1000[i]-1)
}
p1000
p1001

## Probabidade de sair 1010 e 1011

p1010=c(0)
p1011=c(0)
for(i in 1:r){
  p1010[i]= N1010[i] / (N1010[i] + N1011[i])
  p1011[i]=abs(p1010[i]-1)
}
p1010
p1011

## Probabidade de sair 1100 e 1101

p1100=c(0)
p1101=c(0)
for(i in 1:r){
  p1100[i]= N1100[i] / (N1100[i] + N1101[i])
  p1101[i]=abs(p1100[i]-1)
}
p1100
p1101

## Probabidade de sair 1110 e 1111

p1110=c(0)
p1111=c(0)
for(i in 1:r){
  p1110[i]= N1110[i] / (N1110[i] + N1111[i])
  p1111[i]=abs(p1110[i]-1)
}
p1110
p1111

##### Verossimilhança k=3 ####

ver3=c(0)
for(i in 1:r){
  ver3[i]= log2(p0000[i]^N0000[i] * p0001[i]^N0001[i] * p0010[i]^N0010[i] * p0011[i]^N0011[i] * p0100[i]^N0100[i] * p0101[i]^N0101[i] * p0110[i]^N0110[i] * p0111[i]^N0111[i] * p1000[i]^N1000[i] * p1001[i]^N1001[i] * p1010[i]^N1010[i] * p1011[i]^N1011[i] * p1100[i]^N1100[i] * p1101[i]^N1101[i] * p1110[i]^N1110[i] * p1111[i]^N1111[i])
}
#16 MULTIPLICAÇÕES
ver3

##### BIC k=3 ####

bic3=c(0)
for(i in 1:r){
  bic3[i]= ver3[i] - 1/2 * (2^3 * (abs(2)-1) * log(n,2))
}
bic3


######### COMPARANDO OS BIC com k=0 , k=1 , k=2 e k=3 ######### 

cont0=0
cont1=0
cont2=0
cont3=0
for(i in 1:r){
  if(bic0[i]>bic1[i] & bic0[i]>bic2[i] & bic0[i]>bic3[i]){
    cont0=cont0+1
  } else if (bic1[i]>bic0[i] & bic1[i]>bic2[i] & bic1[i]>bic3[i]){
    cont1=cont1+1
  }
  else if(bic2[i]>bic0[i] & bic2[i]>bic1[i] & bic2[i]>bic3[i]){
    cont2 = cont2 +1 
  }
  else{
    cont3=cont3+1
  }
}
cont0
cont1
cont2
cont3

