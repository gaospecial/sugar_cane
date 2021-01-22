
d.all <- fread('data/Compiled Data.txt')
names(d.all) <- gsub(' ', '_', names(d.all))

bootstrp <- function(d) {
    boot <- matrix(0, ncol=6, nrow=5e4)

    d <- as.data.frame(d)
    for (i in 1:5e4) {
        V1.2 <- sample(c(d[,1], d[,2]),replace=T)
        V1.3 <- sample(c(d[,1], d[,3]),replace=T)
        V1.4 <- sample(c(d[,1], d[,4]),replace=T)
        V2.3 <- sample(c(d[,2], d[,3]),replace=T)
        V2.4 <- sample(c(d[,2], d[,4]),replace=T)
        V3.4 <- sample(c(d[,3], d[,4]),replace=T)
        
        boot[i,1] <- mean(V1.2[1:6])-mean(V1.2[7:12])
        boot[i,2] <- mean(V1.3[1:6])-mean(V1.3[7:12])
        boot[i,3] <- mean(V1.4[1:6])-mean(V1.4[7:12])
        boot[i,4] <- mean(V2.3[1:6])-mean(V2.3[7:12])
        boot[i,5] <- mean(V2.4[1:6])-mean(V2.4[7:12])
        boot[i,6] <- mean(V3.4[1:6])-mean(V3.4[7:12])
    }

    aux <- mean(d[,1])-mean(d[,2])
    p1.2 <- sum(aux>boot[,1])/50000

    aux <- mean(d[,1])-mean(d[,3])
    p1.3 <- sum(aux>boot[,2])/50000

    aux <- mean(d[,1])-mean(d[,4])
    p1.4 <- sum(aux>boot[,3])/50000

    aux <- mean(d[,2])-mean(d[,3])
    p2.3 <- sum(aux>boot[,4])/50000

    aux <- mean(d[,2])-mean(d[,4])
    p2.4 <- sum(aux>boot[,5])/50000

    aux <- mean(d[,3])-mean(d[,4])
    p3.4 <- sum(aux>boot[,6])/50000

    ds <- t(combn(4,2))
    r <- data.table(d1 = names(d)[ds[,1]],
                    d1 = names(d)[ds[,2]],
                    pval = c(p1.2, p1.3, p1.4, p2.3, p2.4, p3.4))    
    return(r)
}


ethanol <- bootstrp(d.all[,1:4])
biomass <- bootstrp(d.all[,9:12])
org_ac <- bootstrp(d.all[,13:16])
glycerol <- bootstrp(d.all[,17:20])
viab <- bootstrp(d.all[,21:24])


bootstrp <- function(d) {
    boot <- matrix(0, ncol=6, nrow=5e4)

    d <- as.data.frame(d)
    for (i in 1:5e4) {
        V1.2 <- sample(c(d[-6,1], d[,2]),replace=T)
        V1.3 <- sample(c(d[-6,1], d[,3]),replace=T)
        V1.4 <- sample(c(d[-6,1], d[,4]),replace=T)
        V2.3 <- sample(c(d[,2], d[,3]),replace=T)
        V2.4 <- sample(c(d[,2], d[,4]),replace=T)
        V3.4 <- sample(c(d[,3], d[,4]),replace=T)
        
        boot[i,1] <- mean(V1.2[1:5])-mean(V1.2[6:11])
        boot[i,2] <- mean(V1.3[1:5])-mean(V1.3[6:11])
        boot[i,3] <- mean(V1.4[1:5])-mean(V1.4[6:11])
        boot[i,4] <- mean(V2.3[1:6])-mean(V2.3[7:12])
        boot[i,5] <- mean(V2.4[1:6])-mean(V2.4[7:12])
        boot[i,6] <- mean(V3.4[1:6])-mean(V3.4[7:12])
    }

    aux <- mean(d[-6,1])-mean(d[,2])
    p1.2 <- sum(aux>boot[,1])/50000

    aux <- mean(d[-6,1])-mean(d[,3])
    p1.3 <- sum(aux>boot[,2])/50000

    aux <- mean(d[-6,1])-mean(d[,4])
    p1.4 <- sum(aux>boot[,3])/50000

    aux <- mean(d[,2])-mean(d[,3])
    p2.3 <- sum(aux>boot[,4])/50000

    aux <- mean(d[,2])-mean(d[,4])
    p2.4 <- sum(aux>boot[,5])/50000

    aux <- mean(d[,3])-mean(d[,4])
    p3.4 <- sum(aux>boot[,6])/50000

    ds <- t(combn(4,2))
    r <- data.table(d1 = names(d)[ds[,1]],
                    d1 = names(d)[ds[,2]],
                    pval = c(p1.2, p1.3, p1.4, p2.3, p2.4, p3.4))    
    return(r)
}

prod.all <- bootstrp(d.all[,5:8])

r <- rbind(ethanol, prod.all, biomass, org_ac, glycerol, viab)
fwrite(r, 'all_p_values.csv')
