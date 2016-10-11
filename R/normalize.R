df <- data.frame (one = c(5,2,3,4),
                  two = c(4,1,4,2),
                  three = c(3,4,6,8))

rownames(df) <- toupper(letters[1:4])

df_ranks <- apply(df,2,rank, ties.method = "min")

df_sorted <- data.frame(apply(df,2,sort))
df_up_quart <- apply(df_sorted,1,function(x) quantile(x, prob = .75))

df_normalized <- apply(df_ranks,2, function(x) df_up_quart[x])

rownames(df_normalized) <- rownames(df)

boxplot(df)
boxplot(df_normalized)
