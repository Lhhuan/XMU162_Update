factorial(20)

set.seed(1)
n <- 100
tr <- rbinom(100,1,0.5)
y <- 1 + tr + rnorm(n,0,3)
diff(by(y,tr,mean))
s <- sample(tr,length(tr),FALSE)
diff(by(y,s,mean))

dist <- replicate(2000,diff(by(y,sample(tr,length(tr),FALSE),mean)))






m1 <- lmer(betaval ~ status * cancertype + (1|TSS), control=lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-8)))
  perm <- permanova.lmer(m1,perms = 1E3)
  re <- data.frame(cg = i,status_p = perm$Perm.p[1], cancertype_p = perm$Perm.p[2],
    status_cancertype_p = perm$Perm.p[3])