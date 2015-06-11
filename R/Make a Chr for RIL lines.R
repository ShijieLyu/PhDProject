
pdf("hahaha.pdf")
plot(c(1,10),c(1,10),t="n")
lines(c(2,2),c(4,8),lwd=20,col="red")  # F0, F1
lines(c(3,3),c(4,8),lwd=20,col="blue")  # F0, F1

plot(c(1,15),c(1,10),t="n")    # F2

for(x in 3:10){
  y1 <- runif(1,4,4)
  y2 <- runif(1,5.5,6.8)
  lines(c(x,x),c(4,8),lwd=20,col="red")
  lines(c(x,x),c(y1,y2),lwd=20,col="blue")  
}

plot(c(1,15),c(1,10),t="n")    # F3

for(x in 3:10){
  y1 <- runif(1,4,6)
  y2 <- runif(1,4,6)
  y3 <- runif(1,6,8)
  y4 <- runif(1,6,8)
  lines(c(x,x),c(4,8),lwd=20,col="red")
  lines(c(x,x),c(y1,y2),lwd=20,col="blue") 
  lines(c(x,x),c(y4,y3),lwd=20,col="blue")  
}

#### F0,F1
jpeg("F0.jpg")
plot(c(0,2000),c(0,2500),type = "n")
rect(0,500,100,1900,col="red")
rect(200,500,300,1900,col="blue")
dev.off()
#### F2
plot(c(0,2300),c(0,2500),type = "n")
  for (x in as.numeric(c(0,300,600,900,1200,1500,1800,2100))){
  y1 <- runif(1,500,1200)
  rect(x,500,x+100,y1,col="red")
  rect(x,y1,x+100,1900,col="blue")
}

#### F3

plot(c(0,2300),c(0,2500),type = "n")
  for (x in as.numeric(c(0,300,600,900,1200,1500,1800,2100))){
  y1 <- runif(1,500,1200)
  y2 <- runif(1,600,1900)
  rect(x,500,x+100,1900,col="red")
  rect(x,runif(1,500,700),x+100,runif(1,800,1200),col="blue")
  rect(x,runif(1,1300,1500),x+100,runif(1,1600,1900),col="blue")
}

#### Fn
plot(c(0,2300),c(0,2500),type = "n")
  for (x in as.numeric(c(0,300,600,900,1200,1500,1800,2100))){
  y1 <- runif(1,500,1200)
  rect(x,500,x+100,1900,col="red")
  for(y in 1:40){
    y1 <- runif(1,500,1900)
    rect(x,y1,x+100,y1+25,col="blue")
  }
}



dev.off