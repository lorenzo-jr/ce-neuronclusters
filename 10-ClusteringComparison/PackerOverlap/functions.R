plot.fc <- function(c1, c2)
{
    cp = p[p$cluster == c1,]
    cm = m[m$cluster == c2,]
    gm = as.vector(cm$gene_name)
    gp = as.vector(cp$gene)
    genes = union(gp, gm)

    fc = data.frame(row.names = genes,
        Laurent = rep(0, length(genes)),
        Miller = rep(0, length(genes)))

    fc[gp, 'Laurent'] = cp$avg_logFC
    fc[gm, 'Miller'] = cm$avg_logFC

    plot(fc$Laurent, fc$Miller, xlab = sprintf('Log2 FC cluter %s', c1),
     ylab = sprintf('Log 2 FC cluster %s', c2),
        bg = 'gray', col = 'black',
        pch = 23, xlim = c(-0, 6), ylim = c(-0, 6),
        main = sprintf('%s vs %s', c1, c2))

    gene.sel = genes[order(fc$Laurent, decreasing = T)[1:5]]
    gene.sel = unique(c(gene.sel, genes[order(fc$Miller, decreasing = T)[1:5]]))

    fc.s = fc[gene.sel,]

    points(fc.s$Laurent, fc.s$Miller,
        bg = 'gray', col = 'purple',
        pch = 23)

    text(fc.s$Laurent + 0.3, fc.s$Miller, labels = rownames(fc.s), cex = 0.6, col = 'purple')
}


overlap.pvalue <- function(a, b, c, N = 28000)
{
    m <- a
    n <- N - m
    k <- b
    x <- c
    # return(dhyper(x, m, n, k) + phyper(x, m, n, k, lower.tail = F))
    return(phyper(x, m, n, k, lower.tail = F))
}

intersect.value <- function(a, b, set_size = 26000)
{
    a1 <- length(a)
    a2 <- length(b)
    i  <- length(intersect(a, b))
    c(a = a1, b = a2, i = i, pv = overlap.pvalue(a1, a2, i, set_size))
}