DT <- fread("data/test_data_long")

# remove NA
DT <- na.omit(DT)

# compute price index variable
DT[, price_index := actual_price / regular_price]

# keep only brand A & B & ..
#DT <- DT[brand %chin% c("A","B")] #c("A","B","C", "D", "P")

# take the log
DT[, `:=` (actual_price_log = log(actual_price), price_index_log = log(price_index), volume_sales_log = log(volume_sales))]

# reduce the number of overlapping time series
#DT[brand == "A" & week > 100, `:=` (actual_price_log = NA, price_index_log = NA, volume_sales_log = NA)] #65
#DT[brand == "B" & week < 10, `:=` (actual_price_log = NA, price_index_log = NA, volume_sales_log = NA)] #65

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
X <- as.matrix(DT_input_SUR[, .SD, .SDcols = grep("brand_|actual_price_log_|price_index_log_|FO_|DP_|FD_", names(DT_input_SUR))])

# estimate SUR
m <- itersur(X=X, Y=y, index=index)

# results
m



# next:
#1. use function
#2. improve function

