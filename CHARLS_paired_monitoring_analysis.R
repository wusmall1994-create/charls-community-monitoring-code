# Reproduce the survey-weighted analysis of paired community monitoring in CHARLS.
# Requires R 4.6.1, haven and survey (analysis used survey 4.5).
# Usage: Rscript CHARLS_paired_monitoring_analysis.R <CHARLS root> [output directory]
# Provide an extracted CHARLS directory with the wave folders used below.
# Results are aggregate only. The default output directory is temporary.

if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'Chinese_China.utf8'))
library(haven)
library(survey)
args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript CHARLS_paired_monitoring_analysis.R <CHARLS root> [aggregate output directory]")
b <- args[1]
w <- file.path(b,"2021-2023(WAVE6)")
files <- list.files(w,recursive=TRUE,full.names=TRUE)
read <- function(pattern,cols) {
 f <- grep(pattern,files,value=TRUE)
 stopifnot(length(f)==1)
 read_dta(f,col_select=all_of(cols))
}
s <- read("/Sample_Infor[.]dta$",c("ID","res_wave","died","crosssection"))
wt <- read("/Weights[.]dta$",c("ID","INDV_weight_ad2","INDV_weight_post2"))
v <- c("ID","communityID","xrage","xrgender","dproxy","da011_1_","da012_1_","da011_3_","da012_3_","da035","da136","da036","da137","da013_1_","da014_1_")
d <- do.call(rbind,lapply(2021:2023,function(y){
 x <- read(paste0("/CHARLS-W6-",y,"/Health_Status_and_Functioning[.]dta$"),v)
 x$year <- y
 x
}))
d <- merge(d,s,by="ID")
d <- subset(d,year==res_wave & died==0)
stopifnot(!anyDuplicated(d$ID),!anyDuplicated(wt$ID))
d <- merge(d,wt,by="ID",all.x=TRUE)
p <- read_dta(file.path(b,"2013/CHARLS2013_Dataset/PSU.dta"))
print(table(as_factor(p$urban_nbs)))
stopifnot(!anyDuplicated(p$communityID))
d <- merge(d,p[c("communityID","urban_nbs")],by="communityID",all.x=TRUE)
d$area <- as.character(as_factor(d$urban_nbs))
ht <- d$da011_1_ %in% 1 | d$da012_1_ %in% 1:3
dm <- d$da011_3_ %in% 1 | d$da012_3_ %in% 1:3
d$strict <- ht & dm
d$broad <- (ht | d$da013_1_ %in% 1 | d$da014_1_ %in% 1:3) & dm
d$pair <- d$da035 %in% 1:2 & d$da136 %in% 1:2
d$bp <- as.integer(d$da035==1)
d$diab <- as.integer(d$da136==1)
d$both <- d$bp*d$diab
d$bp_only <- d$bp*(1-d$diab)
d$diab_only <- d$diab*(1-d$bp)
d$neither <- (1-d$bp)*(1-d$diab)
d$gap <- d$bp-d$diab
d$age45 <- d$xrage %in% 45:120
print(c(active=nrow(d),age45=sum(d$age45),age45strict=sum(d$age45&d$strict),age45strictpair=sum(d$age45&d$strict&d$pair)))
print(with(subset(d,age45&strict&pair),table(bp,diab)))
print(with(subset(d,age45&strict&pair),table(crosssection,is.na(INDV_weight_ad2))))
options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=TRUE)
make <- function(wt) {
 x <- d[d$crosssection %in% 1 & is.finite(d[[wt]]) & d[[wt]]>0 & !is.na(d$communityID),]
 x$w <- x[[wt]]
 svydesign(ids=~communityID,weights=~w,data=x)
}
ds <- make("INDV_weight_ad2")
sp <- make("INDV_weight_post2")
main <- subset(ds,age45&strict&pair)
est <- function(z){
 x <- svymean(~bp+diab+both+bp_only+diab_only+neither+gap,z)
 print(c(n=nrow(z$variables),communities=length(unique(z$variables$communityID)),df=degf(z)))
 print(with(z$variables,table(bp,diab)))
 print(round(100*cbind(coef(x),confint(x,df=degf(z))),2))
}
est(main)
for (y in 2021:2023) {cat("YEAR",y,"
"); est(subset(main,year==y))}
for (a in sort(unique(na.omit(main$variables$area)))) {cat("AREA",a,"
"); est(subset(main,area==a))}
print(summary(svyglm(gap~factor(area),subset(main,!is.na(area))))$coefficients)
cat("POST
"); est(subset(sp,age45&strict&pair))
cat("BROAD
"); est(subset(ds,age45&broad&pair))
cat("NO MARKED PROXY
"); est(subset(main,is.na(dproxy)|dproxy!=1))
print(svymean(~xrage+I(xrgender==2),main))
print(table(main$variables$year))
print(table(main$variables$area,useNA="always"))
cat("R version",as.character(getRversion()),"survey",as.character(packageVersion("survey")),"
")
cat("COUNTY SENSITIVITY
")
d$county <- substr(d$communityID,1,4)
print(length(unique(d$county)))
county <- update(ds,county=substr(communityID,1,4))
est(subset(svydesign(ids=~county,weights=~w,data=county$variables),age45&strict&pair))
cat("AREA COMPARISONS
")
print(svymean(~I(area=="Rural"),main))
for(x in c("gap","neither","both")){
 m <- svyglm(reformulate("factor(area)",x),main)
 print(100*c(coef(m)[2],confint(m)[2,]))
}
cat("FREQUENCY
")
print(table(main$variables$da036,useNA="always"))
print(table(main$variables$da137,useNA="always"))

readyear <- function(y,m,cols) read(paste0("/CHARLS-W6-",y,"/",m,"[.]dta$"),cols)
for(y in 2021:2023){
 dem <- readyear(y,"Demographics",c("ID","householdID","xredu"))
 ins <- readyear(y,"Health_Care_Insurance",c("ID","ea001","ea002"))
 hh <- readyear(y,"Household_INC_EXP_Assets",c("householdID","ge003","ge005"))
 health <- readyear(y,"Health_Status_and_Functioning",c("ID","da001","db001","db003","db005","db007","db009"))
 stopifnot(!anyDuplicated(dem$ID),!anyDuplicated(ins$ID),!anyDuplicated(hh$householdID),!anyDuplicated(health$ID))
 a <- merge(merge(merge(dem,ins,by="ID",all.x=TRUE),hh,by="householdID",all.x=TRUE),health,by="ID",all.x=TRUE)
 a$year <- y
 if(y==2021) aux <- a else aux <- rbind(aux,a)
}
dd <- merge(main$variables,aux,by=c("ID","year"),all.x=TRUE)
print(table(as_factor(dd$xredu),useNA="always"))
print(table(dd$ea001,dd$ea002,useNA="always"))
print(table(as_factor(dd$da001),useNA="always"))
print(summary(dd[,c("ge003","ge005")]))
print(with(dd,table(both, freq=da036 %in% 1:7 & da137 %in% 1:7)))
a18 <- file.path(b,"2018/CHARLS2018r")
h18 <- read_dta(file.path(a18,"Health_Status_and_Functioning.dta"),col_select=all_of(c("ID","communityID","da007_1_","da007_3_","zdiagnosed_1_","zdiagnosed_3_","da010_w2_1_1_","da010_w2_2_1_","da010_w2_1_3_","da010_w2_2_3_","da012_w4","da016_w4","da011_w2_1","da011_w4_1")))
demo18 <- read_dta(file.path(a18,"Demographic_Background.dta"),col_select=all_of(c("ID","ba002_1","ba004_w3_1")))
w18 <- read_dta(file.path(a18,"Weights.dta"),col_select=all_of(c("ID","INDV_weight_ad2")))
stopifnot(!anyDuplicated(h18$ID),!anyDuplicated(demo18$ID),!anyDuplicated(w18$ID))
h18 <- merge(merge(h18,demo18,by="ID",all.x=TRUE),w18,by="ID",all.x=TRUE)
h18$age <- 2018-ifelse(h18$ba002_1 %in% 1900:2018,h18$ba002_1,h18$ba004_w3_1)
h18$strict <- (h18$da007_1_ %in% 1 | h18$zdiagnosed_1_ %in% 1) & (h18$da007_3_ %in% 1 | h18$zdiagnosed_3_ %in% 1)
h18$denial <- Reduce("|",lapply(h18[c("da010_w2_1_1_","da010_w2_2_1_","da010_w2_1_3_","da010_w2_2_3_")],function(x)x %in% 99))
h18$strict <- h18$strict & !h18$denial
h18$pair <- h18$da012_w4 %in% 1:2 & h18$da016_w4 %in% 1:2
h18$bp <- as.integer(h18$da012_w4==1)
h18$diab <- as.integer(h18$da016_w4==1)
h18$both <- h18$bp*h18$diab
h18$bp_only <- h18$bp*(1-h18$diab)
h18$diab_only <- h18$diab*(1-h18$bp)
h18$neither <- (1-h18$bp)*(1-h18$diab)
h18$gap <- h18$bp-h18$diab
s18 <- svydesign(ids=~communityID,weights=~INDV_weight_ad2,data=subset(h18,is.finite(INDV_weight_ad2)&INDV_weight_ad2>0))
t18 <- subset(s18,age %in% 45:120 & strict & pair)
cat("2018
")
print(with(subset(h18,age %in% 45:120 & strict),table(pair,useNA="always")))
est(t18)
cat("OVERLAP",sum(t18$variables$ID %in% dd$ID),"
")
dd$area <- factor(dd$area,levels=c("Urban","Rural"))
dd$education <- factor(ifelse(dd$xredu %in% 1:3,"Below primary",ifelse(dd$xredu %in% 4:5,"Primary or middle","High school or above")),levels=c("High school or above","Primary or middle","Below primary"))
dd$insurance <- factor(ifelse(dd$ea001 %in% 3,"None",ifelse(dd$ea002 %in% 1,"Employee",ifelse(dd$ea002 %in% 2,"Resident",NA))),levels=c("Employee","Resident","None"))
dd$food <- ifelse(dd$ge003>0 & dd$ge005>=0,dd$ge005/dd$ge003,NA)
dd$female <- as.integer(dd$xrgender==2)
adl <- as.matrix(dd[c("db001","db003","db005","db007","db009")])
dd$difficulty <- ifelse(rowSums(adl>=2 & adl<=4,na.rm=TRUE)>0,1,ifelse(rowSums(adl==1,na.rm=TRUE)==5,0,NA))

dd$poor <- ifelse(dd$da001 %in% 1:5,as.integer(dd$da001 %in% 4:5),NA)
dd$food100 <- log2(1+dd$food/100)
dd$group <- factor(ifelse(dd$both==1,"Both",ifelse(dd$bp_only==1,"BP only",ifelse(dd$diab_only==1,"Diabetes only","Neither"))),levels=c("Both","BP only","Diabetes only","Neither"))
dd$timing <- factor(ifelse(dd$da036<dd$da137,"BP more often",ifelse(dd$da036>dd$da137,"Diabetes more often","Same")))
ddes <- svydesign(ids=~communityID,weights=~w,data=dd)
print(svymean(~timing,subset(ddes,both==1 & da036 %in% 1:7 & da137 %in% 1:7)))
print(with(subset(dd,both==1),table(timing)))
print(with(subset(h18,age %in% 45:120 & strict),table(is.na(da012_w4),is.na(da016_w4))))
cat("W6 frequency routing
")
print(with(dd,table(bp,da036 %in% 1:7)))
print(with(dd,table(diab,da137 %in% 1:7)))
print(summary(dd$food))
xvars <- c("xrage","female","area","education","insurance","food100","difficulty","poor","year")
print(sapply(dd[xvars],function(x)sum(is.na(x))))
for(g in levels(dd$group)){
 cat("CHAR",g,"
")
 sub <- subset(ddes,group==g)
 print(round(coef(svymean(~xrage+female+area+education+insurance+difficulty+poor,sub,na.rm=TRUE)),3))
 print(svyquantile(~food,sub,quantiles=.5,na.rm=TRUE,ci=FALSE))
}
form <- neither~area+insurance+education+food100+splines::ns(xrage,3)+female+difficulty+poor+factor(year)
fit <- svyglm(form,subset(ddes,complete.cases(dd[xvars])),family=quasibinomial())
print(nobs(fit))
print(table(model.frame(fit)$insurance))
target <- dd[complete.cases(dd[xvars]),]
avg <- function(field,value){
 z <- target
 z[[field]] <- if(is.factor(z[[field]])) factor(value,levels=levels(z[[field]])) else value
 X <- model.matrix(delete.response(terms(fit)),z,contrasts.arg=fit$contrasts)
 pr <- plogis(drop(X %*% coef(fit)))
 list(p=weighted.mean(pr,z$w),grad=colSums(X*(pr*(1-pr)*z$w))/sum(z$w))
}
compare <- function(field,values){
 ref <- avg(field,values[1])
 for(value in values){
  a <- avg(field,value); v <- a$grad-ref$grad
  se <- sqrt(drop(t(v)%*% vcov(fit)%*%v))
  cat(field,value,round(100*c(a$p,a$p-ref$p,(a$p-ref$p)+c(-1,1)*qt(.975,df.residual(fit))*se),2),"
")
 }
}
compare("area",levels(dd$area))
compare("insurance",levels(dd$insurance))
compare("education",levels(dd$education))
compare("food100",log2(1+c(50,150)/100))
out <- if(length(args)>1) args[2] else tempdir()
dir.create(out,showWarnings=FALSE,recursive=TRUE)
line <- function(x) paste0("| ",paste(x,collapse=" | ")," |")
tablemd <- function(x) c(line(colnames(x)),line(rep("---",ncol(x))),apply(x,1,line))
put <- function(x,n) writeLines(enc2utf8(x),file.path(out,n),useBytes=TRUE)
rowstat <- function(z) {
 m <- svymean(~both+bp_only+diab_only+neither+gap,z)
 ci <- confint(m,df=degf(z))*100
 c(nrow(z$variables),sprintf("%.1f (%.1f to %.1f)",100*coef(m),ci[,1],ci[,2]))
}
sets <- list("Wave 6 main"=main,"Age-sex post-stratified"=subset(sp,age45&strict&pair),"Broader hypertension definition"=subset(ds,age45&broad&pair),"Exclude marked proxy responses"=subset(main,is.na(dproxy)|dproxy!=1),"Broader geographic clustering"=subset(svydesign(ids=~county,weights=~w,data=county$variables),age45&strict&pair),"2021 interviews"=subset(main,year==2021),"2022 interviews"=subset(main,year==2022),"2023 interviews"=subset(main,year==2023),"2018 earlier-wave comparison"=t18)
st <- cbind(names(sets),t(vapply(sets,rowstat,character(6))))
colnames(st) <- c("Analysis","n","Both","BP only","Diabetes only","Neither","Difference")
put(tablemd(st),"sensitivity.md")
cf <- function(z) {
 m <- coef(svymean(~xrage+female+area+education+insurance+difficulty+poor,z,na.rm=TRUE))
 c(sprintf("%.1f",m[1]),sprintf("%.1f",100*m[c(2,4:12)]),sprintf("%.1f",coef(svyquantile(~food,z,.5,ci=FALSE,na.rm=TRUE))))
}
ct <- cbind(c("Age, mean (years)","Women (%)","Rural sampling community (%)","High school or above (%)","Primary or middle school (%)","Below primary school (%)","Employee insurance (%)","Resident insurance (%)","No social insurance (%)","Any basic activity difficulty (%)","Poor or very poor health (%)","Food spending per household diner, median (yuan/week)"),sapply(levels(dd$group),function(g)cf(subset(ddes,group==g))))
colnames(ct) <- c("Characteristic","Both (412)","BP only (371)","Diabetes only (59)","Neither (762)")
put(tablemd(ct),"characteristics.md")
ac <- function(field,vals,labels) {
 ref<-avg(field,vals[1])
 t(vapply(seq_along(vals),function(i){
 a<-avg(field,vals[i]); v<-a$grad-ref$grad; V<-vcov(fit); q<-qt(.975,df.residual(fit))
 se<-sqrt(drop(t(v)%*%V%*%v)); sp<-sqrt(drop(t(a$grad)%*%V%*%a$grad))
 fmt<-function(x,s) sprintf("%.1f (%.1f to %.1f)",100*x,100*(x-q*s),100*(x+q*s))
 c(labels[i],fmt(a$p,sp),if(i==1) "Reference" else fmt(a$p-ref$p,se))
 },character(3)))
}
at <- rbind(ac("area",levels(dd$area),c("Urban sampling community","Rural sampling community")),ac("insurance",levels(dd$insurance),c("Employee insurance","Resident insurance","No social insurance")),ac("education",levels(dd$education),levels(dd$education)),ac("food100",log2(1+c(50,150)/100),c("Food spending 50 yuan/week","Food spending 150 yuan/week")))
colnames(at)<-c("Group or value","Adjusted probability, % (95% CI)","Difference, pp (95% CI)")
put(tablemd(at),"adjusted.md")
freq <- subset(ddes,both==1 & da036 %in% 1:7 & da137 %in% 1:7)
fm <- svymean(~timing,freq)
print(round(100*cbind(coef(fm),confint(fm,df=degf(freq))),2))
print(sapply(dd[xvars],function(x)sum(is.na(x))))


denial <- Reduce("|",lapply(h18[c("da010_w2_1_1_","da010_w2_2_1_","da010_w2_1_3_","da010_w2_2_3_")],function(x)x %in% 99))
print(table(denial[h18$strict & h18$age %in% 45:120],h18$pair[h18$strict & h18$age %in% 45:120]))
print(table(is.na(main$variables$INDV_weight_post2)))

print(with(target,table(group)))
z <- subset(t18,!ID %in% dd$ID)
cat("2018 excluding shared respondents
")
est(z)

fmt <- function(x) sprintf("%.1f",x)
desc <- function(z){
 x<-z$variables
 m<-function(v) weighted.mean(x[[v]],x$w,na.rm=TRUE)
 q<-coef(svyquantile(~food,z,c(.25,.5,.75),ci=FALSE,na.rm=TRUE))
 c(fmt(m("xrage")),fmt(100*m("female")),fmt(100*weighted.mean(x$area=="Rural",x$w)),sapply(levels(x$education),function(k)fmt(100*weighted.mean(x$education==k,x$w,na.rm=TRUE))),sapply(levels(x$insurance),function(k)fmt(100*weighted.mean(x$insurance==k,x$w,na.rm=TRUE))),fmt(100*m("difficulty")),fmt(100*m("poor")),sprintf("%.1f (%.1f–%.1f)",q[2],q[1],q[3]))
}
cats<-c(list(ddes),lapply(levels(dd$group),function(g)subset(ddes,group==g)))
ct<-cbind(ct[,1],sapply(cats,desc))
ct[12,1]<-"Food spending, median (IQR), yuan per diner/week"
colnames(ct)<-c("Characteristic","Overall (1604)","Both (412)","BP only (371)","Diabetes only (59)","Neither (762)")
put(tablemd(ct),"characteristics2.md")
stopifnot(abs(sum(sapply(cats[-1],function(z)sum(z$variables$w*z$variables$xrage)))/sum(dd$w)-weighted.mean(dd$xrage,dd$w))<1e-10)
print(sapply(cats,function(z)colSums(is.na(z$variables[c("food","insurance","difficulty","poor")]))))
print(c(zero_diners=sum(dd$ge003==0,na.rm=TRUE),missing_diners=sum(is.na(dd$ge003)),missing_expenditure=sum(is.na(dd$ge005)),negative_expenditure=sum(dd$ge005<0,na.rm=TRUE)))
cc <- complete.cases(dd[xvars])
more <- complete.cases(dd[setdiff(xvars,"poor")])
crud <- function(z,f){
 m<-svyglm(reformulate(f,"neither"),subset(z,!is.na(z$variables[[f]])))
 round(100*cbind(coef(m),confint(m)),2)
}
for(z in list(ddes,subset(ddes,cc))) for(f in c("area","insurance")) print(crud(z,f))
for(keep in list(cc,more)){
 target <- dd[keep,]
 fit <- svyglm(update(form,.~.-poor),subset(ddes,keep),family=quasibinomial())
 cat("Reduced",nobs(fit),"
")
 print(ac("area",levels(dd$area),levels(dd$area)))
 print(ac("insurance",levels(dd$insurance),levels(dd$insurance)))
}
u<-subset(h18,age %in% 45:120 & strict)
missing_bp<-!u$da012_w4 %in% 1:2
missing_dm<-!u$da016_w4 %in% 1:2
skip<-u$da011_w2_1 %in% 2 | u$da011_w4_1 %in% 2
print(table(missing_bp,missing_dm,skip))
for(i in list(u$pair,!u$pair,u$pair & is.finite(u$INDV_weight_ad2)&u$INDV_weight_ad2>0,u$pair & !(is.finite(u$INDV_weight_ad2)&u$INDV_weight_ad2>0))) print(c(n=sum(i),mean_age=mean(u$age[i])))
print(with(dd,tapply(food,year,summary)))
print(with(dd,tapply(ge005,year,function(x)c(n=sum(!is.na(x)),zero=sum(x==0,na.rm=TRUE),max=max(x,na.rm=TRUE),max_count=sum(x==max(x,na.rm=TRUE),na.rm=TRUE)))))
print(with(dd,tapply(ge003,year,summary)))

# Quarterly threshold within the original paired survey domain
stopifnot(with(main$variables, all(bp == 0 | da036 %in% 1:7),
                               all(diab == 0 | da137 %in% 1:7)))
quarterly <- update(main, bp = as.integer(bp == 1 & da036 %in% 1:5),
                         diab = as.integer(diab == 1 & da137 %in% 1:5))
quarterly <- update(quarterly, both = bp * diab, bp_only = bp * (1 - diab),
                    diab_only = diab * (1 - bp), neither = (1 - bp) * (1 - diab),
                    gap = bp - diab)
est(quarterly)
estimate_rows <- function(design, formula) {
  result <- svymean(formula, design)
  100 * cbind(estimate = coef(result), confint(result, df = degf(design)))
}
conditional <- estimate_rows(subset(main, bp == 1), ~I(1-diab))
combined <- estimate_rows(main, ~I(1-neither)+I(bp_only+diab_only))
quarterly_rows <- estimate_rows(quarterly, ~both+bp_only+diab_only+neither+gap)
write.csv(rbind(quarterly_rows, conditional, combined),file.path(out,"quarterly_results.csv"))
