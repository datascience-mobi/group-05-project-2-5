################## variance 

### fold change data

# creat variance data
var_v_FC=as.data.frame(apply(FC,1,var))
var_v_FC

#sort the data
#var_v_FC_min= min(var_v_FC)
var_v_FC_max= max(var_v_FC)

#largest10_var_v_FC_min <- as.data.frame(sort(t(var_v_FC), decreasing = F)[1:10])
largest10_var_v_FC_max <- as.data.frame(sort(t(var_v_FC), decreasing = T)[1:10])

### treated data

# creat variance data
var_TreatedVorinostat=as.data.frame(apply(TreatedVorinostat,1,var))
var_TreatedVorinostat

#sort the data
var_TreatedVorinostat_min= min(var_TreatedVorinostat)
var_TreatedVorinostat_max= max(var_TreatedVorinostat)

largest10_var_TreatedVorinostat_min <- as.data.frame(sort(t(var_TreatedVorinostat), decreasing = F)[1:10])
largest10_var_TreatedVorinostat_max <- as.data.frame(sort(t(var_TreatedVorinostat), decreasing = T)[1:10])
