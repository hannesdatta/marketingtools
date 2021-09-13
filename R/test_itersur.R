library(marketingtools)
library(dplyr)

df <- read.csv('data/test_data_perfect.csv')

# melt data
melted_df <- melt(df, id.vars=c('week','quarter'))
melted_df$brand <- sapply(as.character(melted_df$variable), function(x) rev(strsplit(x, '_', fixed=T)[[1]])[1])
melted_df$measure <- sapply(as.character(melted_df$variable), function(x) paste0(rev(rev(strsplit(x, '_', fixed=T)[[1]])[-1]),collapse='_'))
melted_df$value <- as.numeric(melted_df$value)

df_converted = dcast(melted_df, brand + week ~ measure, measure = 'value', fill = 0)


# add brand fixed effects
for (i in unique(df_converted$brand)) eval(parse(text=paste0('df_converted$dummy_', i,' <- ifelse(df_converted$brand=="', i,'", 1,0)')))


index = df_converted[, c('week','brand')]


X = as.matrix(df_converted[,grep('price|dummy', colnames(df_converted), value=T)])

# get (flattened) Y prices
sales = melted_df[grepl('volume_sales', melted_df$measure),]

y_tmp = merge(index, sales, by = c('brand','week'),all.x=T, sort = FALSE)

y <- as.matrix(as.numeric(y_tmp$value))


m <- itersur(X=X, Y=y, index=index)


