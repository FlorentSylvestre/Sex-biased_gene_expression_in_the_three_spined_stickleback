cov = read.table("coverage.txt")
parse_cov = tapply(cov$V2, cov$V1, sum)

Res = data.frame(Ind = names(parse_cov), Cov = unname(parse_cov))

write.table(x = Res, file = "coverage_sample.txt", row.names =F, col.names = F, quote = F)
