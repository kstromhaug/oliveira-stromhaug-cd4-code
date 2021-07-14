library(tidyverse)



GroupGDClonotypeFamilies <- function(allsplit, gd_clones) {
  all.tcrgroups = data.frame()
  rest = allsplit
  rest$found <-FALSE
  gd_clones$found <- FALSE
  clono.type = 1
  small = 0
  match_rules = c(all=1, insertions=0, deletions=0, substitutions=1)
  show(nrow(gd_clones))
  cont = TRUE
  
  for (i in 1:nrow(gd_clones)) {
    if (i %% 500 == 0){
      cat('row', i, '\n')
    }
    
    cl = gd_clones[i,]$gd_clone
    
    if (cont == TRUE) {
      
      if ( gd_clones[i,]$found==FALSE ) {
        ## PULL OUT ANY IDENTICAL MATCHES, REMOVE FROM LARGER TABLE
        rest[rest$gd_clone==cl, ]$found <- TRUE
        exact.match = rest[rest$gd_clone==cl, ]
        
        if (nrow(exact.match)>0){
          v.alpha <- exact.match[1,]$TRGV_1
          j.alpha <- exact.match[1,]$TRGJ_1
          
          v.alpha.2 <- exact.match[1,]$TRGV_2
          j.alpha.2 <- exact.match[1,]$TRGJ_2
          
          v.beta <- exact.match[1,]$TRDV_1
          d.beta <- exact.match[1,]$TRDD_1
          j.beta <- exact.match[1,]$TRDJ_1
          
          v.beta.2 <- exact.match[1,]$TRDV_2
          d.beta.2 <- exact.match[1,]$TRDD_2
          j.beta.2 <- exact.match[1,]$TRDJ_2
          
          cdr3.a <- exact.match[1,]$TRG_CDR3_1
          cdr3.a.2 <- exact.match[1,]$TRG_CDR3_2
          cdr3.b <- exact.match[1,]$TRD_CDR3_1
          cdr3.b.2 <- exact.match[1,]$TRD_CDR3_2
          
          ## LOOK THROUGH THE REST OF THE gd_cloneS TO FIND ANY IN THE SAME FAMILY
          ## gd_category 1
          if (exact.match[1,]$gd_category=='1.1A1B') {
            ## gd_category 1: match exactly
            vdj.match = rest %>% subset(gd_category=='1.1A1B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) & agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 2: either alpha matches and the beta matches
            vdj.match = rest %>% subset(gd_category=='2.2A1B' & ((TRGV_1==v.alpha & TRGJ_1==j.alpha) | (TRGV_2==v.alpha & TRGJ_2==j.alpha)) & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset((agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) | agrepl(cdr3.a, TRG_CDR3_2, max.distance = match_rules)) & agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 3: alpha matches and one of betas matches
            vdj.match = rest %>% subset(gd_category=='3.1A2B' & (TRGV_1==v.alpha & TRGJ_1==j.alpha) & ((TRDV_2==v.beta & TRDJ_2==j.beta) | (TRDV_1==v.beta & TRDJ_1==j.beta)) & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) & (agrepl(cdr3.b, TRD_CDR3_2, max.distance = match_rules) | agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules)))
              exact.match = rbind(exact.match, cdr3.match)
            }
            ## gd_category 4: match alpha but no beta
            vdj.match = rest %>% subset(gd_category=='4.A0B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 5: match beta but no alpha
            vdj.match = rest %>% subset(gd_category=='5.B0A' & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 6: two alphas match but no beta
            vdj.match = rest %>% subset(gd_category=='6.2A0B' & ((TRGV_1==v.alpha & TRGJ_1==j.alpha) | (TRGV_2==v.alpha & TRGJ_2==j.alpha)) & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) | agrepl(cdr3.a, TRG_CDR3_2, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 7: two betas match but no alpha
            vdj.match = rest %>% subset(gd_category=='7.2B0A' & ((TRDV_1==v.beta & TRDJ_1==j.beta) | (TRDV_2==v.beta & TRDJ_2==j.beta)) & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules) | agrepl(cdr3.b, TRD_CDR3_2, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            rest[rest$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            gd_clones[gd_clones$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            
          } else if (exact.match[1,]$gd_category=='2.2A1B') {
            ## gd_category 1: match exactly
            vdj.match = rest %>% subset(gd_category=='2.2A1B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & TRGV_2==v.alpha.2 & TRGJ_2==j.alpha.2 & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) & agrepl(cdr3.a.2, TRG_CDR3_2, max.distance = match_rules) & agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 2: first alpha and beta matches
            vdj.match = rest %>% subset(gd_category=='1.1A1B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) & agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 3: second alpha and beta matches
            vdj.match = rest %>% subset(gd_category=='1.1A1B' & TRGV_1==v.alpha.2 & TRGJ_1==j.alpha.2 & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a.2, TRG_CDR3_1, max.distance = match_rules) & agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 4: an alpha matches but no beta
            vdj.match = rest %>% subset(gd_category=='4.A0B' & ((TRGV_1==v.alpha & TRGJ_1==j.alpha) | (TRGV_1==v.alpha.2 & TRGJ_1==j.alpha.2)) & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) | agrepl(cdr3.a.2, TRG_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 6: two alphas match but no beta
            vdj.match = rest %>% subset(gd_category=='6.2A0B' & ((TRGV_1==v.alpha & TRGJ_1==j.alpha) & (TRGV_2==v.alpha.2 & TRGJ_2==j.alpha.2)) & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) & agrepl(cdr3.a.2, TRG_CDR3_2, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## gd_category 5: match beta but no alpha
            vdj.match = rest %>% subset(gd_category=='5.B0A' & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
              # rest = rest[!(rest$gd_clone %in% cdr3.match$gd_clone), ]
            }
            
            # show(nrow(rest[rest$gd_clone %in% clonotype.group$gd_clone, ]))
            rest[rest$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            gd_clones[gd_clones$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            # clonotype.group$clonotype <- clono.type
            
          } else if (exact.match[1,]$gd_category=='3.1A2B') {
            ## gd_category 1: match exactly
            vdj.match = rest %>% subset(gd_category=='3.1A2B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & TRDV_2==v.beta.2 & TRDJ_2==j.beta.2 & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) & agrepl(cdr3.b.2, TRD_CDR3_2, max.distance = match_rules) & agrepl(cdr3.b, TRD_CDR3_1, max.distance = c(1,1,1,1)))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## gd_category 2: first alpha and beta matches
            vdj.match = rest %>% subset(gd_category=='1.1A1B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) & agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## gd_category 3: first alpha and second beta matches
            vdj.match = rest %>% subset(gd_category=='1.1A1B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & TRDV_1==v.beta.2 & TRDJ_1==j.beta.2 & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules) & agrepl(cdr3.b.2, TRD_CDR3_1, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## gd_category 4: alpha matches but no betas
            vdj.match = rest %>% subset(gd_category=='4.A0B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## gd_category 5: either beta matches but no alpha
            vdj.match = rest %>% subset(gd_category=='5.B0A' & ((TRDV_1==v.beta & TRDJ_1==j.beta) | (TRDV_2==v.beta.2 & TRDJ_2==j.beta.2)) & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules) | agrepl(cdr3.b.2, TRD_CDR3_2, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## gd_category 7: two betas match but no alpha
            vdj.match = rest %>% subset(gd_category=='7.2B0A' & ((TRDV_1==v.beta & TRDJ_1==j.beta) & (TRDV_2==v.beta.2 & TRDJ_2==j.beta.2)) & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules) & agrepl(cdr3.b.2, TRD_CDR3_2, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            rest[rest$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            gd_clones[gd_clones$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            
          } else if (exact.match[1,]$gd_category=='4.A0B') {
            ## gd_category 4: alpha matches but no betas
            vdj.match = rest %>% subset(gd_category=='4.A0B' & TRGV_1==v.alpha & TRGJ_1==j.alpha & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, TRG_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            rest[rest$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            gd_clones[gd_clones$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            
          } else if (exact.match[1,]$gd_category=='5.B0A') {
            ## gd_category 4: alpha matches but no betas
            vdj.match = rest %>% subset(gd_category=='5.B0A' & TRDV_1==v.beta & TRDJ_1==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, TRD_CDR3_1, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            rest[rest$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
            gd_clones[gd_clones$gd_clone %in% exact.match$gd_clone, ]$found<-TRUE
          } # end gd_category checking
        } 
        
        if (nrow(exact.match) > 0) {
          exact.match$gd_clonotype <- clono.type
          clono.type = clono.type + 1
        } 
        all.tcrgroups = rbind(all.tcrgroups, exact.match)
        gd_clones[i,]$found<-TRUE
      }
    }
  }
  
  return(all.tcrgroups)
}
