library(data.table)
library(systemfit)

DT <- fread("data/test_data_long.csv")

# remove NA
DT <- na.omit(DT)

# compute price index variable
DT[, price_index := actual_price / regular_price]

# keep only brand A & B & ..
DT <- DT[brand %chin% c("A","B", "C", "P")] #c("A","B","C", "D", "P")

# take the log
DT[, `:=` (actual_price_log = log(actual_price), price_index_log = log(price_index), volume_sales_log = log(volume_sales))]

# reduce the number of overlapping time series
#DT[brand == "A" & week > 65, `:=` (actual_price_log = NA, price_index_log = NA, volume_sales_log = NA)] #65 #64
#DT[brand == "B" & week < 65, `:=` (actual_price_log = NA, price_index_log = NA, volume_sales_log = NA)] #65
# delete NA values
#DT <- na.omit(DT)

# remove redundant columns
DT[, c("regular_price", "quarter", "FD_FO_DO", "actual_price", "price_index", "volume_sales") := NULL]

# stacked X data.table
DT_input_SUR <- data.table::dcast(DT, brand + week + volume_sales_log ~ brand, value.var = c("actual_price_log", "price_index_log", "FO", "DO", "FD"), fill = 0)
# add dummies
inds <- unique(DT_input_SUR$brand)
DT_input_SUR[, paste0("brand_", inds) := lapply(inds, function(x) as.numeric(brand == x))]

# create the objects: (i) index (ii) y and (iii) X. These are needed as input for the function itersur.
index <- DT_input_SUR[, .(week, brand)]
y <- as.matrix(DT_input_SUR[, volume_sales_log])
X <- as.matrix(DT_input_SUR[, .SD, .SDcols = grep("brand_|actual_price_log_|price_index_log_|FO_|DO_|FD_", names(DT_input_SUR))]) #|FO_|DO_|FD_

# estimate SUR
m <- itersur(X=X, Y=y, index=index, method="FGLS-Praise-Winsten") #method="FGLS-Praise-Winsten"

# results
m

## compare with systemfit
#DT_wide <- dcast(DT, week ~ brand, value.var = c("volume_sales_log", "actual_price_log", "price_index_log", "FO", "DO", "FD"))
#EqA <- volume_sales_log_A ~ actual_price_log_A + price_index_log_A + FO_A + DO_A + FD_A
#EqB <- volume_sales_log_B ~ actual_price_log_B + price_index_log_B + FO_B + DO_B + FD_B
#EqC <- volume_sales_log_C ~ actual_price_log_C + price_index_log_C + FO_C + DO_C + FD_C
#EqP <- volume_sales_log_P ~ actual_price_log_P + price_index_log_P + FO_P + DO_P + FD_P

#system <- list(EqA = EqA, EqB = EqB, EqC = EqC, EqP = EqP)
#SUR_systemfit <- systemfit(system, data = DT_wide, method = "SUR", maxiter = 500)
#summary(SUR_systemfit)




# next:
#1. use function
#2. improve function

# function need to be robust to:
# 1. NA values in the data
# 2. few variation in the data
# 3. no/little overlap in time-series

