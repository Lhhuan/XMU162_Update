interval<-seq(from=500000, to=1000000, by=500000)

for (j in interval){
    dir_name = paste("int",j, sep = "")
    if(dir.exists(dir_name)){
        print(dir_name exit)
    }
    else{
        dir.create(dir_name)
    }
}