                            
                            # TOPIC 1
#read data
a = read.table("67.txt")

#get the lines with sequence only
linesWithSeq = grep(pattern="^>", x = a[,1], invert = TRUE)

#strsplit each sequence to handle them easier
data = strsplit(a[linesWithSeq,],"")

      # QUESTION 1 

usefuldata = list()

counter = 1

#get usefuldata (all sequences with 6 or more positions )
for (i in 1:length(data)){
  
  if (length(data[[i]]) >= 6){
    
    usefuldata[[counter]] = data[[i]]
    counter = counter + 1
  }
} 


length(data)-length(usefuldata) # only 2 sequences with less than 6 elements


substrings = list()
index = list()  
len = list ()
idx = 1

#function to get all the substrings with length 6 that have haming distance <2 from the motif.
Hammingdistance = function(x){
  
  # Given motif
  motif = c("C","G","T","C","A","C")
  
  dif = vector(mode = "numeric") 
  
  # for each sequence substrings with length 6 will be in total:
  #length(sequence) - length(motif) + 1
  
  for (i in 1:(length(x) - length(motif) + 1)){
    
    # x[i:(i+5)] is the current substring
    
    #hamming distance between the motif  and the substring
    dif = sum(motif != x[i:(i+5)])
    
    #if hamming distance <2:
    if (dif == 0 | dif == 1){
      
      # get the substring
      substrings[[idx]] <<- x[i:(i+5)] 
      
      # get substrings starting position
      index[[idx]] <<- i 
      
      # get the length of the sequence  that this substring belongs to
      len[[idx]] <<- length(x)
      
      idx <<- idx + 1
    }
  }
  
}

# apply that function to usefuldata 
ham=lapply(usefuldata,Hammingdistance) 


      # QUESTION 2 : Construct the PWM using these substrings

mtrx = t(matrix(unlist(substrings) , nrow = 6))

pfm = matrix(0, nrow=4, ncol=6)

rownames(pfm) = c("A","C","G","T")

# count the number of appearances   of a letter in each position
pfm = apply(mtrx, 2, function(x){table(factor(x, levels= c("A","C","G","T")))})

# make count to frequency
ppm = pfm/as.vector(apply(pfm, 2, sum))

#  divide by 1/4 and take the log2!
pwm = log2(ppm/0.25)

pwm

# QUESTION 3 : Construct the logo for the PPM

require(ggplot2)
require(ggseqlogo)

ggseqlogo(ppm,method='bits')
Sys.sleep(1)

#ggseqlogo(pwm,ncol=4)


        # QUESTION 4 : 
#How far is the starting  position of each substring  from the end of the sequence

position_from_theend = unlist(len) - unlist(index)

plot(density(position_from_theend))



# What do you observe?

# see the distribution for the length of the sequences
plot(density(unlist(len)))

# percentage of sequences with less than 1000 positions
length( which(unlist(len)<1000) )/length(unlist(len))

#Only 1.1% of the sequences have length less than 1000.

#That means that, if distances from the end were random, 
#the distribution of the distances would follow the uniform distribution in range [5,1000]

#Instead we see  a right skew in the distribution of the distances,
#meaning that there are more substrings (with hamming distance < 2 from the motif)
# starting positions closer to the end of the sequences.



                # TOPIC 2

#read lines of our dataset
dataset = readLines("GDS3713.soft")

#remove all lines that start with ^, #, !
cleanlines = dataset[grep('^[/!#/^]', dataset, invert = TRUE)]

#make new file with our clean dataset
write(cleanlines,"GDS3713.soft.clean",ncolumns = length(cleanlines[1]), sep="\t")

data=read.table("GDS3713.soft.clean",header = TRUE)

dim(data)

#remove the first 2 columns that are not expression values !!!
data=data[,-c(1,2)]

#log10 the data
data=log10(data)


# function to apply on every row(gene) of data to get the p value

#according to the experiment info, first 40 samples are control and the next 39 are smokers!
getpval<- function(row){
  
  control_i=row[1:40]
  
  smokers_i=row[41:79]
  
  # get pvalue for ??1: smokers > control  (t test)
  pvalue = t.test(smokers_i,control_i,alternative = 'greater')$p.value
  return(pvalue)
}

#get the pvalues
pvalues = apply(data, 1, getpval)

# apply fdr correction
pvalues.fdr=p.adjust(pvalues,method = "fdr")

## genes with adjusted pvalues<0.05 have value statistically higher in smokers than non-smokers

## get the genes with adjusted pvalues<0.05 !!
genes_higher_smokers=which(pvalues.fdr<0.05)










                        #TOPIC 3

# Read ms files!!
read.ms.output <- function( file.ms.output=NA ) {
  
  txt=NA  
  
  if( !is.na(file.ms.output) ) txt <- scan(file=file.ms.output,
                                           what="character", sep="\n", quiet=TRUE)
  if( is.na(txt[1]) ){
    print("Usage: read.ms.output(txt), or read.ms.output(file=filename)")
    return()
  }
  nsam <- as.integer( strsplit(txt[1], split=" ")[[1]][2] )
  ndraws <- as.integer( strsplit( txt[1], split=" ")[[1]][3] )
  
  h <- numeric()
  result <- list()
  gamlist <- list()
  positions <- list()
  
  marker <- grep("prob",txt)
  probs <- sapply(strsplit(txt[marker], split=":"), function(vec) as.numeric(vec[2]))
  marker <- grep("time",txt)
  times <- sapply(strsplit(txt[marker], split="\t"), function(vec){ as.numeric(vec[2:3])} )
  
  
  ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
  marker <- grep("segsites", txt)
  
  if( length(marker) != ndraws){
    stop( paste("length: ", length(marker), " ndraws: ", ndraws) )
    stopifnot(length(marker) == ndraws)
  }
  
  ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
  segsites <- sapply(strsplit(txt[marker], split=" "), function(vec) as.integer(vec[2]) )
  for(draw in seq(along=marker)) {
    if(!(draw %% 100)) cat(draw, " ")
    if(segsites[draw] > 0) {
      tpos <- strsplit(txt[marker[draw]+1], split=" ")
      positions[[draw]] <- as.numeric( tpos[[1]][ 2:(segsites[draw]+1) ] ) 
      haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
      haplotypes <- strsplit(haplotypes, split="")
      h <- sapply(haplotypes, function(el) c(as.integer(el)))
      ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
      if(segsites[draw] == 1) h <- as.matrix(h)
      ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
      else h <- t(h)
    }
    else {
      h <- matrix(nrow=nsam, ncol=0)
      positions[[draw]]<- NA
    }
    gamlist[[draw]] <- h
    stopifnot(all(dim(h) == c(nsam, segsites[draw]))) 
  }
  cat("\n")
  list(segsites=segsites, gametes=gamlist, probs=probs, times=t(times), positions=positions, nsam=nsam, nreps=ndraws ) 
}





#   Get simulation data

#download.file("http://139.91.162.101//teaching/project2022/R/ms/ms.sim.out", "ms.sim.out")
simdata<- read.ms.output("ms.sim.out")



# CASE 1  (2 summary statistics: Spns , mean pairwise differences )

reps=length(simdata$gametes)

#matrix for the summary statistics
all.sim.stats = matrix(0, nrow=reps, ncol=2) 


# function to calculate number of SNPs (Single-nucleotide polymorphism)
snp=function(data){
  
  snps=0
  for (i in 1:(ncol(data))){
    #an uparxei 1 sthn sthlh tote einai polumorfikh
    if(length(which(data[,i]==1))>0){
      snps=snps+1
    }
  }
  return(snps)
}

#save SNPs as the first summary statistic
all.sim.stats[,1]=sapply(simdata$gametes,snp)


#function to calculate the average number of pairwise differences between sequences
pairdif = function(data){
  dif = 0
  for (i in 1:(nrow(data)-1)){
    for(j in (i+1):(nrow(data))){
      dif = dif + sum(data[i,] != data[j,])
    }
  }
  return(2*dif/(nrow(data)*(nrow(data)-1)))
}

#save pairwise differences as the second summary statistic
all.sim.stats[,2]=as.numeric(lapply(simdata$gametes, pairdif)) 


## for the observation

#download.file("http://139.91.162.101//teaching/project2022/R/ms/ms.obs.out","ms.obs.out")
#read the observation data:
ms <- read.ms.output("ms.obs.out")
data = ms$gametes[[1]]

#statistics of the observation 
stats = c(snp(ms$gametes[[1]]), pairdif(ms$gametes[[1]])) 

#euclidian distance
eucl = function(m, obs){
  difs = apply(m, 1, function(x){ sqrt(sum((obs - x)^2)) })
  return(difs)
}

#euclidian distances for the statics of the observation and the simulations
difs = eucl(all.sim.stats, stats)

#closest 500 simulations to observation
accepted.dif.indexes  = which(difs <= sort(difs)[500] ) 

#download.file("http://139.91.162.101//teaching/project2022/R/ms/params.txt","params.txt")


# params=G values of the simulations
params=read.table("params.txt")

# G values for the closest 500 simulations
posterior=params[accepted.dif.indexes,]

#prior density
d.prior=density(params[,1])

#posterior density
d.posterior=density(posterior)

# our estimate for the growth rate G
case1_infer=mean(posterior)

#plot prior and posterior 
plot(d.prior, col='1',ylim=c(0,0.03),main="Case 1", xlab="Growth Rate")
points(d.posterior, col='6', type='l')
legend('topleft',legend=c("prior",'posterior'),fill=c(1,6))


library('abc')
myabc = abc( target = stats , param = params , sumstat = all.sim.stats , method = 'loclinear' ,tol= 0.0501)
abc1=myabc$unadj.values


points(density(abc1),type='l')
mean(abc1)
case1_infer


# CASE 2 with Tajima's D statistic

all.sim.stats2 = matrix(0, nrow=reps, ncol=3)

#first 2 statistics are the same
all.sim.stats2[,1:2]=all.sim.stats

# for Tajimas D:

#D=(average_pairwise_diferences- Spns/normalization)/normalization2

#D=(khat- S/a1)/normalization2

#from wikipedia:
tajimasd=function (x,khat,S){
  n=nrow(x)
  tmp <- 1:(n - 1)
  a1 <- sum(1/tmp)
  a2 <- sum(1/tmp^2)
  b1 <- (n + 1)/(3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  D <- (khat - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
  return(D)
}


#save tajimas D as the 3rd statistic
all.sim.stats2[,3]=mapply(tajimasd,simdata$gametes,khat=all.sim.stats2[,2] ,S=all.sim.stats2[,1])


# stats for the observation
stats2= c(snp(ms$gametes[[1]]), pairdif(ms$gametes[[1]]), tajimasd(ms$gametes[[1]],stats[2],stats[1])) 

# euclidian distances for the statics of the observation and the simulations
difs2 = eucl(all.sim.stats2, stats2)

# closest 500 simulations to observation
accepted.dif.indexes2= which(difs2 <= sort(difs2)[500] ) 

# G values for the closest 500 simulations
posterior2=params[accepted.dif.indexes2,]

#posterior density
d.posterior2=density(posterior2)

# our estimate for the growth rate G
case2_infer=mean(posterior2)
case2_infer

#plot prior and posterior 
plot(d.prior, col='1',ylim=c(0,0.03),main="Case 2", xlab="G")
points(d.posterior2, col='6', type='l')
legend('topleft',legend=c("prior",'posterior'),fill=c(1,6))

# 
# library('abc')
# myabc2 = abc( target = stats2 , param = params , sumstat = all.sim.stats2 , method = 'loclinear' ,tol= 0.0501)
# abc2=myabc2$unadj.values
# 
# 
# points(density(abc2),type='l')
# 
# mean(abc2)
# case2_infer
# 


