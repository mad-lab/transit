from rpy2.robjects import r, globalenv

def def_r_clustering():
    r('''
    r_clustering = function(data_filename, K, D) {
        # input: temp_LFCs.txt
        data = read.table(data_filename,head=T,sep='\t')
        lfcs = data[,3:length(colnames(data))]
        labels = paste(data$Rv,'/',data$Gene,sep="")
        rownames(lfcs) = labels

        R = length(rownames(lfcs)) # genes
        C = length(colnames(lfcs)) # conditions

        library(corrplot)

        fname = "lfcs_boxplot.png"
        png(fname,width=300+15*C,height=500)
        boxplot(lfcs,las=2,ylab="log2 fold change of genes (relative to the median)")
        dev.off()

        km = kmeans(lfcs,K)
        clusters = km$cluster

        write.table(data.frame(lfcs,clust=clusters),"temp_clust.txt", sep='\t', quote=F)

        pca = prcomp(t(lfcs))

        library(MASS)
        mat = as.matrix(lfcs)
        ldapca = lda(clusters~.,lfcs)
        X = mat %*% ldapca$scaling[,1]
        Y = mat %*% ldapca$scaling[,2]
        fname = "pca_genes.png"
        png(fname,width=1000,height=1000)
        plot(X,Y,pch=20)
        text(X,Y,label=labels,adj=c(0.5,1.5),col=clusters)
        dev.off()

        actpca = prcomp(lfcs,center=TRUE,scale=TRUE)
        fname = "pca_conditions.png"
        png(fname,width=500,height=500)
        plot(actpca$rotation[,1:2],pch=20)#,xlim=c(-1,1),ylim=c(-1,1))
        text(actpca$rotation[,1:2],label=colnames(lfcs),adj=c(0.5,1.5),cex.lab=1.5)
        dev.off()

        print(length(lfcs[,1]))
        hc = hclust(dist(lfcs),method="ward.D2")
        print("****hclust_genes****")
        print(hc)
        fname = "hclust_genes.png"
        png(fname,width=300+15*R,height=600)
        plot(hc)
        K2 = max(2,min(K,max(hc$height)))
        rect.hclust(hc,k=K2,border='red')
        dev.off()

        hc = hclust(dist(t(lfcs)),method="ward.D2")
        print("****hclust_conditions****")
        print(hc)
        fname = "hclust_conditions.png"
        png(fname,width=300+15*C,height=600)
        # plot(hc)
        # D2 = max(2,min(D,C-1))
        # rect.hclust(hc,k=D2,border='red')
        # dev.off()

        fname = "cluster_opt.png"
        png(fname)
        library(factoextra)
        fviz_nbclust(lfcs, kmeans, method = "wss",k.max=30) # or silhouette or gap_stat
        #fviz_nbclust ( lfcs, kmeans, method = "wss",k.max=30)+geom_hline ( yintercept=52 )+geom_vline(xintercept=10)
        dev.off()

        ####################################################
        # required factoextra and corrplot libraries

        var = get_pca_var(actpca) 
        fname = "condition_PCs.png"
        png(fname,width=1000,height=1000)
        corrplot(var$cor,main="Principle Components",mar=c(1,1,1,1))
        dev.off()

        S = diag(actpca$sdev,D,D)
        rawLoadings = actpca$rotation[,1:D] %*% S
        vmax = varimax(rawLoadings)
        rotatedLoadings = vmax$loadings
        rotatedScores = scale(actpca$x[,1:D]) %*% vmax$rotmat

        # https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r
        # scores = U.S = X.V = actpca$x = scale(lfcs) %*% actpca$rotation
        combined_transform = (actpca$rotation[,1:D] %*% S) %*% vmax$rotmat
        # vmax$loadings is same as combined_transform, but with abs(vals)<0.1 blanked out
        write.table(round(combined_transform,6),"varimax_loadings.txt",sep='\t',quote=F)

        fname = "varimax.png"
        png(fname,width=800,height=800)
        corrplot(rotatedLoadings,main="Varimax loadings",mar=c(1,1,1,1)) 
        dev.off()

        scores = rotatedScores
        colnames(scores) = paste("score",seq(1:D),sep="")
        squares = scores*scores
        ss = apply(squares,1,sum)
        cos2 = squares/ss
        best = as.matrix(apply(cos2,1,max))
        assoc = cos2*sign(scores)
        colnames(assoc) = paste("assoc",seq(1:D),sep="")

        # # generate null distribution for cos2 to closest axis
        # 
        # library(MASS)
        # SAMPLES = 10000
        # sample = mvrnorm(n=SAMPLES,mu=rep(0,D),Sigma=diag(D))
        # # no need to apply Varimax rotation
        # squares = sample*sample
        # ss = apply(squares,1,sum)
        # Xcos2 = squares/ss
        # Xbest = as.matrix(apply(Xcos2,1,max))
        # #distn = ecdf(best)
        # pvals = apply(best,1,function (x) { length(Xbest[Xbest>=x]) })
        # pvals = pvals/SAMPLES
        # padj = p.adjust(pvals,method="BH")
         
        #res = data.frame(data,scores,assoc,best,pvals,padj)
        res = data.frame(data,scores,assoc)
        write.table(format(res,digits=3,scientific=F), "temp_scores.txt",sep='\t',quote=F,row.names=F)

    }
    ''')
    return globalenv['r_clustering']

def def_r_samples_corrplot():
    r('''
    r_samples_corrplot = function(data_filename) {
        data = read.table(data_filename, sep='\t',head=T)
        vals = as.matrix(data[,4:length(colnames(data))])
        N = length(colnames(vals))

        library(corrplot)

        fname = "samples_corrplot.png"
        png(fname,width=300+20*N,height=300+20*N)
        #corrplot(cor(vals)) # among TAsites (not good)

        TAsites = aggregate ( coord~ORF,data=data,length)
        colnames(TAsites) = c("ORF","sites")
        temp = aggregate(vals,by=list(data$ORF),FUN=mean)
        corrplot(cor(temp[TAsites$sites>=10,2:length(colnames(temp))]))

        dev.off()
    }''')

    return globalenv['r_samples_corrplot']

def def_r_conditions_corrplot():
    r('''
        r_conditions_corrplot = function(data_filename) {
            data = read.table(data_filename, sep='\t', head=T)
            vals = as.matrix(data[,3:length(colnames(data))])
            N = length(colnames(vals))

            library(corrplot)

            png("conditions_corrplot.png", width=20*N+300, height=20*N+300)
            corrplot(cor(vals))
            dev.off()
        }
    ''')

    return globalenv['r_conditions_corrplot']


def def_r_make_heatmap():
    r('''
        r_make_heatmap = function(data_filename) {
            data = read.table(data_filename, head=T, sep='\t')
            lfcs = as.matrix(data[,3:length(colnames(data))])
            labels = paste(data$Rv,'/',data$Gene,sep="")
            rownames(lfcs) = labels

            library(corrplot)
            library(RColorBrewer)
            library(gplots)

            # Red for down regulated, green for upregulated. 
            redgreen <- function(n) { c( hsv(h=0/6, v=seq(1,0,length=n/2) ), hsv(h=2/6, v=seq(0,1,length=n/2) ), ) }
            colors <- colorRampPalette(c("red", "white", "green"))(n = 1000)

            C = length(colnames(lfcs))
            R = length(rownames(lfcs))
            W = 300+C*30
            H = 300+R*15

            fname = "heatmap.png"
            png(fname,width=W,height=H)
            #par(oma=c(10,1,1,10)) # b,l,t,r; units="lines of space"
            heatmap.2(lfcs,col=colors,margin=c(12,12),lwid=c(1,7),lhei=c(1,7),trace="none",cexCol=1.4,cexRow=1.4) # make sure white=0
            dev.off()
        }
    ''')
    return globalenv['r_make_heatmap']

