#### data.table notes
### try use column names rather then numbers/index to access columns
### instead of using DT[,1:5], use DT[,columnname1:columnname5]
### or DT[,c("columnname1","columnname2")]
### more preferred DT[,.(columnname1,columnname2)],use column names directly as variables within .()
### you can place any expression, oprations within .()
### DT[1,] or DT[,1] will return a data.table object rather than vector
### select a single column by name DT$columnname or DT[["columnname"]] will return a vector
### DT[,columnname] also return a vector, DT[,.(columnname)] will return a data.table
### DT[i,j,by]

dataMat <- data.table::fread('GroupAK.somaticvariants.tsv',sep = '\t',header = T,check.names = F)

## rename columns
setnames(dataMat,'#CHROM','CHROM')

## get first 5 columns
dataMat[,1:5]
dataMat[,CHROM:ALT]
dataMat[,.(CHROM,POS,ID,REF,ALT)]

## get column based on value from a variable
tmp <- c('CHROM','POS')
dataMat[,..tmp]
dataMat[,tmp,with=FALSE]

## DT[i,j,by], you can write any expression within i or j
## DT[i=raw,j=column,by=groupby]
dataMat[CHROM=='X',.(newc=paste0(CHROM,'_',POS))]
## .SD is keyword meaning rest of columns excluding the grouping variable
dataMat[ , lapply(.SD, unique), by = CHROM]
## for a subset of .SD, you can specify .SDcols
dataMat[ , lapply(.SD, unique), by = CHROM,.SDcols=4:5]
dataMat[ , sum(POS), by = CHROM]

## create a new data.table using another as template
dataMat[0]

## perform join between two data.tables
X = data.table(grp = c("a", "a", "b",
                       "b", "b", "c", "c"), foo = 1:7)
setkey(X, grp)
Y = data.table(grp=c("b","b","c", "c"), bar = c(4,3,2,1))
setkey(Y,grp)
X[Y,mult='all']
X[Y,mult='first']
##INNER join
X[Y, nomatch=0]
merge(X, Y, all=FALSE)
##LEFT OUTER	
Y[X]	
merge(X, Y, all.x=TRUE)
##RIGHT OUTER	
X[Y]	merge(X, Y, all.y=TRUE)
## FULL OUTER
merge(X, Y, all=TRUE)


## special bulid-in variable .N
## number of observations in current group
flights[, .(.N), by = .(origin)]
## sort by ascending or descending order
flights[order(origin, -dest)]

## return output sorted
flights[carrier == "AA",
        .(mean(arr_delay), mean(dep_delay)),
        keyby = .(origin, dest, month)]



