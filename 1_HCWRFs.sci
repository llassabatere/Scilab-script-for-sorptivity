// Description
    
    tic();  // clock
    
    // Computation of sorptivity and c_p parameter by L. Lassabatere et al. (2021)
    /// long computations --> please wait!!

// Initialization
  
    clear
    mode(0);
    ieee(1);
    clearglobal;
  
  // root directory for Scilab computations
  
    T_wd = pwd();
    
    T_rep = T_wd + "\0_functions";
    
    chdir(T_rep);
    
    exec('0_general_tool_functions.sci',-1);
    exec('1_WRCHCF_func.sci',-1);
    
    close_windows(20);
    
    i_fig = 1;
    
    chdir(T_wd);
    
///////// values of saturation degrees and water potentials ///////////////////

    xh_star = logspace(-6,6,101)';
    xSe = linspace(0,1,201)';
    xh_star_lim = logspace(16,-6,501)';
    
    xR = linspace(0.025,0.975,20)';
    // values of shape parameters
    
    
//////// Synthetic soils !! BC !! functions /////////////////////////////////////////
    
 // values of shape parameters //
    
    lambda_BC = 2*xR./(1-xR);         // in agreement with the values considered for vGB
    eta_BC = 2./lambda_BC+3;
    
    ntot_BC = size(lambda_BC,1);
    nleg_BC = floor(ntot_BC);
    
 // color scale for plots // 
    
    cmap_BC = rainbowcolormap(nleg_BC);
    i0 = 0;
    for i = 1:nleg_BC
        id_color_BC(nleg_BC-i+1) = color(cmap_BC(i+i0,1)*255,cmap_BC(i+i0,2)*255,cmap_BC(i+i0,3)*255);
    end
    
// WR dimensionless functions
    
    ySe_BC = cell();
    yh_star_BC = cell();
    
    for i = 1:ntot_BC
        ySe_BC{i} = Se_hstar_BC(xh_star,lambda_BC(i));                           // Fig subplot(221) - Se(h)
        yh_star_BC{i} = hstar_Se_BC(xSe,lambda_BC(i));
    end
    
    yKr_h_BC = cell();
    yKr_Se_BC = cell();
    yKr_Se_bis_BC = cell();
    xSe_Kr_bis_BC = cell();
    yh_star_Kr_bis_BC = cell();

    for i = 1:ntot_BC
        yKr_Se_BC{i} = Kr_BC(xSe,eta_BC(i));
        
        [a b] = find(yKr_Se_BC{i} == 0);                        // remove zero values
        nnKr_Se_BC(i) = max(a)+1;
        yKr_Se_bis_BC{i} = yKr_Se_BC{i}(nnKr_Se_BC(i):$);
        xSe_Kr_bis_BC{i} = xSe(nnKr_Se_BC(i):$);
        yh_star_Kr_bis_BC{i} = hstar_Se_BC(xSe_Kr_bis_BC{i},lambda_BC(i));
        Se_Kh_lim_BC(i) = xSe(nnKr_Se_BC(i));
        
        [a b] = find(yKr_Se_BC{i} > 10^-15);                  // remove small values
        yKr_Se_bis2_BC{i} = yKr_Se_BC{i}(a);
        xSe_Kr_bis2_BC{i} = xSe(a);
        yh_star_Kr_bis2_BC{i} = hstar_Se_BC(xSe_Kr_bis2_BC{i},lambda_BC(i));
        Se_Kh_lim2_BC(i) = xSe_Kr_bis2_BC{i}(1);
        
    end
   
    for i = 1:ntot_BC
        yKr_h_BC{i} = Kr_BC(Se_hstar_BC(xh_star_lim,lambda_BC(i)),eta_BC(i));    // Fig subplot(222) - K(Se)
        
        [a b] = find(yKr_h_BC{i} == 0);
        if a == [] then a = 0;end;
        nnKr_h_BC(i) = max(a)+1;
        yKr_h_bis_BC{i} = yKr_h_BC{i}(nnKr_h_BC(i):$);
        xh_star_Kr_bis_BC{i} = xh_star_lim(nnKr_h_BC(i):$);
        xh_Kh_lim_BC(i) = xh_star_lim(nnKr_h_BC(i));
        
        [a b] = find(yKr_h_BC{i} > 10^-15);
        yKr_h_bis2_BC{i} = yKr_h_BC{i}(a);
        xh_star_Kr_bis2_BC{i} = xh_star_lim(a);
        ySe_Kr_bis2_BC{i} = Se_hstar_BC(xh_star_Kr_bis2_BC{i},lambda_BC(i));      // Fig subplot(222) - K(h)
    end

    yDstar_h_BC = cell();
    yDstar_Se_BC = cell();
    yDstar_Se_bis_BC = cell();
    xSe_Dstar_bis_BC = cell();
    yh_star_Dstar_bis_BC = cell();

    for i = 1:ntot_BC
        yDstar_Se_BC{i} = abs(D_etoile_BC(xSe,lambda_BC(i),eta_BC(i)));           // Fig subplot(224) - D(Se)
        
        [a b] = find(yDstar_Se_BC{i} == 0);
        if a == [] then a = 1;end;
        nnDstar_BC(i) = max(a)+1;
        yDstar_Se_bis_BC{i} = yDstar_Se_BC{i}(nnDstar_BC(i):$);
        xSe_Dstar_bis_BC{i} = xSe(nnDstar_BC(i):$);
        yh_star_Dstar_bis_BC{i} = yh_star_BC{i}(nnDstar_BC(i):$);
        xSe_Dstar_lim_BC(i) = xSe(nnDstar_BC(i));
        
        [a b] = find(yDstar_Se_BC{i} > 10^-15);
        yDstar_Se_bis2_BC{i} = yDstar_Se_BC{i}(a);
        xSe_Dstar_bis2_BC{i} = xSe(a);
        yh_star_Dstar_bis2_BC{i} = yh_star_BC{i}(a);
    end
    
//////// Synthetic soils !! vGB !! functions /////////////////////////////////////////
    
 // values of shape parameters //
    
    m_vGB = xR;
    n_vGB = floor(2./(1-m_vGB)*100)/100;
    eta_vGB = 2./(n_vGB-2)+3;
    
    ntot_vGB = size(n_vGB,1);
    nleg_vGB = floor(ntot_vGB);
    
 // color scale for plots // 
    
    cmap_vGB = rainbowcolormap(nleg_vGB);
    i0 = 0;
    for i = 1:nleg_vGB
        id_color_vGB(nleg_vGB-i+1) = color(cmap_vGB(i+i0,1)*255,cmap_vGB(i+i0,2)*255,cmap_vGB(i+i0,3)*255);
    end
    
// WR dimensionless functions
    
    ySe_vGB = cell();
    yh_star_vGB = cell();
    
    for i = 1:ntot_vGB
        ySe_vGB{i} = Se_hstar_vGB(xh_star,n_vGB(i));                             // Fig subplot(221) - Se(h)
        yh_star_vGB{i} = hstar_Se_vGB(xSe,n_vGB(i));
    end

    yKr_h_vGB = cell();
    yKr_Se_vGB = cell();
    yKr_Se_bis_vGB = cell();
    xSe_Kr_bis_vGB = cell();
    yh_star_Kr_bis_vGB = cell();

    for i = 1:ntot_vGB
        yKr_Se_vGB{i} = Kr_vGB(xSe,eta_vGB(i));
        
        [a b] = find(yKr_Se_vGB{i} == 0);
        nnKr_Se_vGB(i) = max(a)+1;
        yKr_Se_bis_vGB{i} = yKr_Se_vGB{i}(nnKr_Se_vGB(i):$);
        xSe_Kr_bis_vGB{i} = xSe(nnKr_Se_vGB(i):$);
        yh_star_Kr_bis_vGB{i} = hstar_Se_vGB(xSe_Kr_bis_vGB{i},n_vGB(i));
        Se_Kh_lim_vGB(i) = xSe(nnKr_Se_vGB(i));
        
        [a b] = find(yKr_Se_vGB{i} < 10^-16);
        nnKr_Se2_vGB(i) = max(a)+1;
        yKr_Se_bis2_vGB{i} = yKr_Se_vGB{i}(nnKr_Se2_vGB(i):$);
        xSe_Kr_bis2_vGB{i} = xSe(nnKr_Se2_vGB(i):$);
        yh_star_Kr_bis2_vGB{i} = hstar_Se_vGB(xSe_Kr_bis2_vGB{i},n_vGB(i));
        Se_Kh_lim2_vGB(i) = xSe(nnKr_Se2_vGB(i));
        
    end
   
    for i = 1:ntot_vGB
        yKr_h_vGB{i} = Kr_vGB(Se_hstar_vGB(xh_star_lim,n_vGB(i)),eta_vGB(i));     // Fig subplot(223) - K(h)
        
        [a b] = find(yKr_h_vGB{i} == 0);
        if a == [] then a = 0;end;
        nnKr_h_vGB(i) = max(a)+1;
        yKr_h_bis_vGB{i} = yKr_h_vGB{i}(nnKr_h_vGB(i):$);
        xh_star_Kr_bis_vGB{i} = xh_star_lim(nnKr_h_vGB(i):$);
        xh_Kh_lim_vGB(i) = xh_star_lim(nnKr_h_vGB(i));
        
        [a b] = find(yKr_h_vGB{i} > 10^-10);
        yKr_h_bis2_vGB{i} = yKr_h_vGB{i}(a);
        xh_star_Kr_bis2_vGB{i} = xh_star_lim(a);
        ySe_Kr_bis2_vGB{i} = Se_hstar_vGB(xh_star_Kr_bis2_vGB{i},n_vGB(i));      // Fig subplot(223) - K(Se)
        
    end

    yDstar_h_vGB = cell();
    yDstar_Se_vGB = cell();
    yDstar_Se_bis_vGB = cell();
    xSe_Dstar_bis_vGB = cell();
    yh_star_Dstar_bis_vGB = cell();

    for i = 1:ntot_vGB
        yDstar_Se_vGB{i} = abs(D_etoile_vGB(xSe,n_vGB(i),eta_vGB(i)));           // Fig subplot(223) - D(Se)
        
        [a b] = find(yDstar_Se_vGB{i} == 0);
        if a == [] then a = 1;end;
        nnDstar_vGB(i) = max(a)+1;
        yDstar_Se_bis_vGB{i} = yDstar_Se_vGB{i}(nnDstar_vGB(i):$);
        xSe_Dstar_bis_vGB{i} = xSe(nnDstar_vGB(i):$);
        yh_star_Dstar_bis_vGB{i} = yh_star_vGB{i}(nnDstar_vGB(i):$);
        xSe_Dstar_lim_vGB(i) = xSe(nnDstar_vGB(i));
        
        [a b] = find(yDstar_Se_vGB{i} > 10^-10);
        yDstar_Se_bis2_vGB{i} = yDstar_Se_vGB{i}(a);
        xSe_Dstar_bis2_vGB{i} = xSe(a);
        yh_star_Dstar_bis2_vGB{i} = yh_star_vGB{i}(a);
    end
    
//////// Synthetic soils !! vGB !! functions /////////////////////////////////////////
    
 // values of shape parameters //
    
    m_vGM = xR;
    n_vGM = floor(1./(1-m_vGM)*100)/100;
    
    l_vGM = 0.5*ones(n_vGM);
    
    ntot_vGM = size(n_vGM,1);
    nleg_vGM = floor(ntot_vGM);
    
 // color scale for plots // 
    
    cmap_vGM = rainbowcolormap(nleg_vGM);
    i0 = 0;
    for i = 1:nleg_vGM
        id_color_vGM(nleg_vGM-i+1) = color(cmap_vGM(i+i0,1)*255,cmap_vGM(i+i0,2)*255,cmap_vGM(i+i0,3)*255);
    end
    
// WR dimensionless functions
    
    ySe_vGM = cell();
    yh_star_vGM = cell();
    
    for i = 1:ntot_vGM
        ySe_vGM_vGM{i} = Se_hstar_vGM(xh_star,n_vGM(i));                          // Fig subplot(221) - Se(h)
        yh_star_vGM{i} = hstar_Se_vGM(xSe,n_vGM(i));
    end
    
    yKr_h_vGM = cell();
    yKr_Se_vGM = cell();
    yKr_Se_bis_vGM = cell();
    xSe_Kr_bis_vGM = cell();
    yh_star_Kr_bis_vGM = cell();

    for i = 1:ntot_vGM
        yKr_Se_vGM{i} = Kr_vGM(xSe,n_vGM(i),l_vGM(i));                           // Fig subplot(223) - K(Se)
        
        [a b] = find(yKr_Se_vGM{i} == 0);
        nnKr_Se_vGM(i) = max(a)+1;
        yKr_Se_bis_vGM{i} = yKr_Se_vGM{i}(nnKr_Se_vGM(i):$);
        xSe_Kr_bis_vGM{i} = xSe(nnKr_Se_vGM(i):$);
        yh_star_Kr_bis_vGM{i} = yh_star_vGM{i}(nnKr_Se_vGM(i):$);
        Se_Kh_lim_vGM(i) = xSe(nnKr_Se_vGM(i));
        
        [a b] = find(yKr_Se_vGM{i} < 10^-15);
        nnKr_Se2_vGM(i) = max(a)+1;
        yKr_Se_bis2_vGM{i} = yKr_Se_vGM{i}(nnKr_Se2_vGM(i):$);
        xSe_Kr_bis2_vGM{i} = xSe(nnKr_Se2_vGM(i):$);
        yh_star_Kr_bis2_vGM{i} = yh_star_vGM{i}(nnKr_Se2_vGM(i):$);
        Se_Kh_lim2_vGM(i) = xSe(nnKr_Se2_vGM(i));
        
    end
   
    for i = 1:ntot_vGM
        yKr_h_vGM{i} = Kr_vGM(Se_hstar_vGM(xh_star_lim,n_vGM(i)),n_vGM(i),l_vGM(i));     // Fig subplot(223) - K(h)
        
        [a b] = find(yKr_h_vGM{i} == 0);
        if a == [] then a = 0;end;
        nnKr_h_vGM(i) = max(a)+1;
        yKr_h_bis_vGM{i} = yKr_h_vGM{i}(nnKr_h_vGM(i):$);
        xh_star_Kr_bis_vGM{i} = xh_star_lim(nnKr_h_vGM(i):$);
        xh_Kh_lim_vGM(i) = xh_star_lim(nnKr_h_vGM(i));
        
        [a b] = find(yKr_h_vGM{i} > 10^-15);
        yKr_h_bis2_vGM{i} = yKr_h_vGM{i}(a);
        xh_star_Kr_bis2_vGM{i} = xh_star_lim(a);
        ySe_Kr_bis2_vGM{i} = Se_hstar_vGM(xh_star_Kr_bis2_vGM{i},n_vGM(i));
    end

    yDstar_h_vGM = cell();
    yDstar_Se_vGM = cell();
    yDstar_Se_bis_vGM = cell();
    xSe_Dstar_bis_vGM = cell();
    yh_star_Dstar_bis_vGM = cell();

    for i = 1:ntot_vGM
        
        yDstar_Se_vGM{i} = abs(D_etoile_vGM(xSe,n_vGM(i),l_vGM(i)));             // Fig subplot(223) - D(Se)
        
        [a b] = find(yDstar_Se_vGM{i} == 0);
        if a == [] then a = 1;end;
        nnDstar_vGM(i) = max(a)+1;
        yDstar_Se_bis_vGM{i} = yDstar_Se_vGM{i}(nnDstar_vGM(i):$);
        xSe_Dstar_bis_vGM{i} = xSe(nnDstar_vGM(i):$);
        yh_star_Dstar_bis_vGM{i} = yh_star_vGM{i}(nnDstar_vGM(i):$);
        xSe_Dstar_lim_vGM(i) = xSe(nnDstar_vGM(i));
        
        [a b] = find(yDstar_Se_vGM{i} > 10^-15);
        yDstar_Se_bis2_vGM{i} = yDstar_Se_vGM{i}(a);
        xSe_Dstar_bis2_vGM{i} = xSe(a);
        yh_star_Dstar_bis2_vGM{i} = yh_star_vGM{i}(a);
    end
    
    
//////// Synthetic soils !! KG !! functions /////////////////////////////////////////
    
 // values of shape parameters //
    
    sigma_KG = 1./xR -1;                 // fix the same values for shape indexes
    l_KG = 0.5*ones(sigma_KG);
    
    ntot_KG = size(sigma_KG,1);
    nleg_KG = floor(ntot_KG);
    
 // color scale for plots // 
    
    cmap_KG = rainbowcolormap(nleg_KG);
    i0 = 0;
    for i = 1:nleg_KG
        id_color_KG(i) = color(cmap_KG(i+i0,1)*255,cmap_KG(i+i0,2)*255,cmap_KG(i+i0,3)*255);
    end
    
// WR dimensionless functions
    
    nn_KG_min = 3;                         // to avoid boundaries - small or very hight values of m
    nn_KG_max = ntot_KG-2;
    
    ySe_KG = cell();
    yh_star_KG = cell();
    
    for i = nn_KG_min:nn_KG_max
        ySe_KG{i} = Se_hstar_KG(xh_star,sigma_KG(ntot_KG-i+1));                  // Fig subplot(221) - Se(h)
        yh_star_KG{i} = hstar_Se_KG(xSe,sigma_KG(ntot_KG-i+1));
    end
    
    
    yKr_h_KG = cell();
    yKr_Se_KG = cell();
    yKr_Se_bis_KG = cell();
    xSe_Kr_bis_KG = cell();
    yh_star_Kr_bis_KG = cell();

    for i = nn_KG_min:nn_KG_max
        yKr_Se_KG{i} = Kr_KG(xSe,sigma_KG(ntot_KG-i+1),l_KG(ntot_KG-i+1));       // Fig subplot(222) - K(Se)
        
        [a b] = find(yKr_Se_KG{i} == 0);
        nnKr_Se_KG(i) = max(a)+1;
        yKr_Se_bis_KG{i} = yKr_Se_KG{i}(nnKr_Se_KG(i):$);
        xSe_Kr_bis_KG{i} = xSe(nnKr_Se_KG(i):$);
        yh_star_Kr_bis_KG{i} = yh_star_KG{i}(nnKr_Se_KG(i):$);
        Se_Kh_lim_KG(i) = xSe(nnKr_Se_KG(i));
        
        [a b] = find(yKr_Se_KG{i} > 10^-15);
        nnKr_Se2_KG(i) = max(a)+1;
        yKr_Se_bis2_KG{i} = yKr_Se_KG{i}(a);
        xSe_Kr_bis2_KG{i} = xSe(a);
        yh_star_Kr_bis2_KG{i} = hstar_Se_KG(xSe_Kr_bis2_KG{i},sigma_KG(i)); 
        Se_Kh_lim2_KG(i) = xSe_Kr_bis2_KG{i}(1);
    end
   
    for i = nn_KG_min:nn_KG_max
        ySe_2 = Se_hstar_KG(xh_star_lim,sigma_KG(ntot_KG-i+1));
        
        yKr_h_KG{i} = Kr_KG(ySe_2,sigma_KG(ntot_KG-i+1),l_KG(ntot_KG-i+1));      // Fig subplot(223) - K(h)
        
        [a b] = find(yKr_h_KG{i} == 0);
        if a == [] then a = 0;end;
        nnKr_h_KG(i) = max(a)+1;
        yKr_h_bis_KG{i} = yKr_h_KG{i}(nnKr_h_KG(i):$);
        xh_star_Kr_bis_KG{i} = xh_star_lim(nnKr_h_KG(i):$);
        xh_Kh_lim_KG(i) = xh_star_lim(nnKr_h_KG(i));
        
        [a b] = find(yKr_h_KG{i} > 10^-15);
        yKr_h_bis2_KG{i} = yKr_h_KG{i}(a);
        xh_star_Kr_bis2_KG{i} = xh_star_lim(a);
        ySe_Kr_bis2_KG{i} = ySe_2(a);
    end

    yDstar_h_KG = cell();
    yDstar_Se_KG = cell();
    yDstar_Se_bis_KG = cell();
    xSe_Dstar_bis_KG = cell();
    yh_star_Dstar_bis_KG = cell();

    for i = nn_KG_min:nn_KG_max
        yDstar_Se_KG{i} = abs(D_etoile_KG(xSe,sigma_KG(ntot_KG-i+1),l_KG(ntot_KG-i+1)));    // Fig subplot(224) - D(Se)
        
        [a b] = find(yDstar_Se_KG{i} == 0);
        if a == [] then a = 1;end;
        nnDstar_KG(i) = max(a)+1;
        yDstar_Se_bis_KG{i} = yDstar_Se_KG{i}(nnDstar_KG(i):$);
        xSe_Dstar_bis_KG{i} = xSe(nnDstar_KG(i):$);
        yh_star_Dstar_bis_KG{i} = yh_star_KG{i}(nnDstar_KG(i):$);
        xSe_Dstar_lim_KG(i) = xSe(nnDstar_KG(i));
        
        [a b] = find(yDstar_Se_KG{i} > 10^-15);
        yDstar_Se_bis2_KG{i} = yDstar_Se_KG{i}(a);
        xSe_Dstar_bis2_KG{i} = xSe(a);
        yh_star_Dstar_bis2_KG{i} = yh_star_KG{i}(a);
    end
    
    
//////////////////////////   Figure of the paper   ////////////////// 
    
  // legends
    
    legtxt_BC_ = "$"+string(floor(m_vGB*1000)/1000)+"$";
    legtxt_KG_ = legtxt_BC_;
    legtxt_vGM_ = legtxt_BC_;
    legtxt_vGB_ = legtxt_BC_;
    
    nn_BC = [[1:2:size(lambda_BC,1)]'(1:$-1);size(lambda_BC,1)];
    
    legtxt_BC = legtxt_BC_(nn_BC);

    
///////// Graphical options ///////////////////////////////////////////////////
  
    xf1 = 3;               // font_size axis graduations
    xf2 = 4;                  // font_size axis title
    xleg = 3;
    
    xfs_alphabet = 4;
    
    alphabet = "$\bf {\large {("+["a" "b" "c" "d" "e" "f" "g" ..
                    "h" "i" "j" "k" "l" "m" "n" ..
                    "o" "p" "q" "r" "s" "t" "u" ..
                    "v" "w" "x" "y" "z"]+")}}$";
    
////////////////////////////// margins ////////////////////

 //   amargins = [margin_left,margin_right,margin_top,margin_bottom];
 
     xx_margins = [0.15,0.1,0.1,0.2];
    
////////////////////////////// WRHCFs ////////////////////
    
    xm1 = 10^-3;
    xM1 = 10^7;
    ym1 = 0;
    yM1 = 1;
    
    xm2 = 0;
    xM2 = 1;
    ym2 = 10^-14;
    yM2 = 10^1;
    
    
// //// Fig 1 //////////////////////////////////////
    
    clf(i_fig);
    scf(i_fig);
    
    nfig_WRHCFs = i_fig;
    
    i_fig = i_fig+1;
    
    fig = gcf();
    fig.figure_size = [900,788];
    fig.auto_resize = "on";
    
    legtxt_fig1 = "$\textrm{"+["BC";"vGB";"vGM";"KG"]+"}$";
    
    
 /// Plot Se(h) // 
    
    xnf = [1:4];
    i0 = 6;                                                // selected curve
    
    disp("Selected curves for Figure 1");
    disp("xR = "+string(xR(i0)));
    disp("lambda_BC = "+string(lambda_BC(i0)));
    disp("n_vGB = "+string(n_vGB(i0)));
    disp("n_vGM = "+string(n_vGM(i0)));
    disp("sigma_KG = "+string(sigma_KG(i0)));
    

    ii = 1;
    
    subplot(2,2,xnf(ii))
    
        plot(yh_star_BC{i0}(2:($-1)),xSe(2:($-1)),"blue -");
        plot(yh_star_vGB{i0}(2:($-1)),xSe(2:($-1)),"red -");
        plot(yh_star_vGM{i0}(2:($-1)),xSe(2:($-1)),"color","green");
        plot(yh_star_KG{ntot_BC-i0}(2:($-1)),xSe(2:($-1)),"color","black");
        
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [xm1,ym1;xM1,yM1];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lnn";
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        for i=1:4
            a.children(i).children.thickness = 1.3;
        end
        
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$S_{e}$","fontsize",xf2);
        
        leg = legend(legtxt_fig1,pos=4,boxed=%f);
        leg.font_size = xleg;
        leg.fill_mode = "off";
        
        xstring(xM1,0.9,alphabet(xnf(ii)));
        
        ii = ii + 1;
       
        
 /// Plot K(Se) //
        
    subplot(2,2,xnf(ii))
    
        plot(xSe_Kr_bis_BC{i0},yKr_Se_bis_BC{i0},"blue -");
        plot(xSe_Kr_bis_vGB{i0},yKr_Se_bis_vGB{i0},"red --");
        plot(xSe_Kr_bis_vGM{i0},yKr_Se_bis_vGM{i0},"color","green");
        plot(xSe_Kr_bis_KG{ntot_BC-i0},yKr_Se_bis_KG{ntot_BC-i0},"color","black");
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nln";
        a.data_bounds = [xm2,ym2;xM2,yM2];
        a.tight_limits = ["on","on","off"];
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        for i=1:4
            a.children(i).children.thickness = 1.3;
        end
        
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$K_{r}$","fontsize",xf2);
          
        xstring(0,0.3,alphabet(xnf(ii)));
        
        ii = ii + 1;
        
        
 /// Plot K(h) //

    subplot(2,2,xnf(ii))
    
        plot(xh_star_Kr_bis2_BC{i0},yKr_h_bis2_BC{i0},"blue -");
        plot(xh_star_Kr_bis2_vGB{i0},yKr_h_bis2_vGB{i0},"red -");
        plot(xh_star_Kr_bis2_vGM{i0},yKr_h_bis2_vGM{i0},"color","green");
        plot(xh_star_Kr_bis2_KG{ntot_BC-i0},yKr_h_bis2_KG{ntot_BC-i0},"color","black");
        
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [10^-3,10^-8;10^6,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lln";
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        for i=1:4
            a.children(i).children.thickness = 1.3;
        end
    
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$K_{r}$","fontsize",xf2);
        
        xstring(8*10^5,1,alphabet(xnf(ii)));
        
        ii = ii+1;

 /// plots of D(Se) //


    subplot(2,2,xnf(ii))
    
        plot(xSe_Dstar_bis2_BC{i0}(1:($-1)),yDstar_Se_bis2_BC{i0}(1:($-1)),"blue -");
        plot(xSe_Dstar_bis2_vGB{i0}(1:($-1)),yDstar_Se_bis2_vGB{i0}(1:($-1)),"red -");
        plot(xSe_Dstar_bis2_vGM{i0}(1:($-1)),yDstar_Se_bis2_vGM{i0}(1:($-1)),"color","green");
        plot(xSe_Dstar_bis2_KG{ntot_BC-i0}(1:($-1)),yDstar_Se_bis2_KG{ntot_BC-i0}(1:($-1)),"color","black");
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [0,10^-8;1,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "nln";
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        for i=1:4
            a.children(i).children.thickness = 1.3;
        end
          
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$D^*$","fontsize",xf2);

        xstring(0,1.5,alphabet(xnf(ii)));

        ii = ii+1;
//        
    
 ////// Fig 2 //////////////////////////////////////
//

    clf(i_fig);
    scf(i_fig);
    
    nfig_WRHCFs_xr = i_fig;
    
    i_fig = i_fig+1;
    
    fig = gcf();
    fig.figure_size = [900,788];
    fig.auto_resize = "on";
    
    xx = [1:4:13]';
    xnf = [xx;xx+1;xx+2;xx+3];
    
    ii = 1;
    
// Options
    
    xar = 1.5;
    
    ry_abt = 0.95;
    ry_abt2 = 0.1;
  
    xf1 = 0.25;               // font_size axis graduations
    xf2 = 3;                  // font_size axis title
    xleg = 0.1;
    
    xfs_alphabet = 1.5;
    
    alphabet = "$\textrm{("+["a" "b" "c" "d" "e" "f" "g" ..
                    "h" "i" "j" "k" "l" "m" "n" ..
                    "o" "p" "q" "r" "s" "t" "u" ..
                    "v" "w" "x" "y" "z"]+")}$";
 
 
     xx_margins = [0.25,0.1,0.1,0.27];
     
    
 /// Plot Se(h) /////////////////////////////////////////////
    
    subplot(4,4,xnf(ii))
    
        for i = 1:size(nn_BC,1)
            plot(yh_star_BC{nn_BC(i)}(2:($-1)),xSe(2:($-1)),"color",id_color_BC(nn_BC(i)));
        end
        
        for i = 1:ntot_BC
            plot(yh_star_BC{i}(2:($-1)),xSe(2:($-1)),"color",id_color_BC(i));
        end
        
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$S_{e,\,BC}$","fontsize",xf2);
        
        leg = legend(legtxt_BC,pos=4,boxed=%f);
        leg.font_size = xleg;
        leg.fill_mode = "off";
    
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [xm1,ym1;xM1,yM1];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lnn";
        a.font_size = xf1;
//        atl = a.x_ticks.labels;
//        a.x_ticks.labels = "-"+atl;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
        xarrows([xM1/10;xM1/10],[Se_hstar_BC(xM1/10,lambda_BC(1));Se_hstar_BC(xM1/10,lambda_BC($))],10*xar);
        
        xstring(xM1,0.84*ry_abt,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        
        ii = ii + 1;
      
    subplot(4,4,xnf(ii))
        
        for i = 1:ntot_vGB
            plot(yh_star_vGB{i}(2:($-1)),xSe(2:($-1)),"color",id_color_vGB(i));
        end
        
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$S_{e,\,vGB}$","fontsize",xf2);
        
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [xm1,ym1;xM1,yM1];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lnn";
        a.font_size = xf1;
//        atl = a.x_ticks.labels;
//        a.x_ticks.labels = "-"+atl;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
        xarrows([xM1/10;xM1/10],[Se_hstar_vGB(xM1/10,n_vGB(1));Se_hstar_vGB(xM1/10,n_vGB($))],10*xar);
        
        xstring(xM1,0.84*ry_abt,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        
        ii = ii + 1;
        
    subplot(4,4,xnf(ii))
        
        for i = 1:ntot_vGM
            plot(yh_star_vGM{i}(2:($-1)),xSe(2:($-1)),"color",id_color_vGM(i));
        end
        
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$S_{e,\,vGM}$","fontsize",xf2);
        
        xarrows([xM1/10;xM1/10],[Se_hstar_vGM(xM1/10,n_vGM(1));Se_hstar_vGM(xM1/10,n_vGM($))],10*xar);
        
        xstring(xM1,0.84*ry_abt,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [xm1,ym1;xM1,yM1];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lnn";
        a.font_size = xf1;
//        atl = a.x_ticks.labels;
//        a.x_ticks.labels = "-"+atl;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
        ii = ii + 1;
      
    subplot(4,4,xnf(ii))
        
        for i = nn_KG_min:nn_KG_max
            plot(yh_star_KG{i}(2:($-1)),xSe(2:($-1)),"color",id_color_KG(i));
        end
        
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$S_{e,\,KG}$","fontsize",xf2);
        
        xarrows([xm1*10;xm1*10],[Se_hstar_KG(xm1*10,sigma_KG(3));Se_hstar_KG(xm1*10,sigma_KG($-2))],10*xar);
        
        xstring(xM1,0.84*ry_abt,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [xm1,ym1;xM1,yM1];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lnn";
        a.font_size = xf1;
//        atl = a.x_ticks.labels;
//        a.x_ticks.labels = "-"+atl;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
        ii = ii + 1;
        
 /// Plot K(Se) ///////////////////////////////////////////// 
    
    subplot(4,4,xnf(ii))
    
        for i = 1:ntot_BC
            plot(xSe_Kr_bis_BC{i},yKr_Se_bis_BC{i},"color",id_color_BC(i));
        end
        
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$K_{r,\,BC}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nln";
        a.data_bounds = [xm2,ym2;xM2,yM2];
        a.tight_limits = ["on","on","off"];
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
        xarrows([0.75;0.75],[Kr_BC(0.75,eta_BC(1));Kr_BC(0.75,eta_BC($))],12.5*xar);
          
        xstring(0,yM2/(0.3*10^2)*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        
        ii = ii + 1;
    
    subplot(4,4,xnf(ii))
    
        for i = 1:ntot_vGB
            plot(xSe_Kr_bis_vGB{i},yKr_Se_bis_vGB{i},"color",id_color_vGB(i));
        end
        
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$K_{r,\,vGB}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nln";
        a.data_bounds = [xm2,ym2;xM2,yM2];
        a.tight_limits = ["on","on","off"];
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
        xarrows([0.75;0.75],[Kr_vGB(0.75,eta_vGB(1));Kr_vGB(0.75,eta_vGB($))],12.5*xar);
        
        xstring(0,yM2/(0.3*10^2)*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        
        ii = ii + 1;
    
    subplot(4,4,xnf(ii))
    
        for i = 1:ntot_vGM
            plot(xSe_Kr_bis_vGM{i},yKr_Se_bis_vGM{i},"color",id_color_vGM(i));
        end
          
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$K_{r,\,vGM}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
       
        
        xarrows([0.75;0.75],[Kr_vGM(0.75,n_vGM(1),l_vGM(1));Kr_vGM(0.75,n_vGM($),l_vGM($))],12.5*xar);
        
        xstring(0,yM2/(0.3*10^2)*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nln";
        a.data_bounds = [xm2,ym2;xM2,yM2];
        a.tight_limits = ["on","on","off"];
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
        ii = ii + 1;
    
    subplot(4,4,xnf(ii))
    
        for i = nn_KG_min:nn_KG_max
            plot(xSe_Kr_bis_KG{i},yKr_Se_bis_KG{i},"color",id_color_KG(i));
        end
        
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$K_{r,\,KG}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
    
        xarrows([0.75;0.75],[Kr_KG(0.75,sigma_KG(3),l_KG(3));Kr_KG(0.75,sigma_KG($-2),l_KG($-2))],12.5*xar);
        
        xstring(0,yM2/(0.3*10^2)*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nln";
        a.data_bounds = [xm2,ym2;xM2,yM2];
        a.tight_limits = ["on","on","off"];
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
        ii = ii + 1;
        
        
 /// Plot K(h) ///////////////////////////////////////////// 

// position of the arrows for K(h)
    
    [a1 b1] = min(abs(10^-7-yKr_h_bis2_BC{1}));
    [a2 b2] = min(abs(10^-7-yKr_h_bis2_BC{$}));
    [a3 b3] = min(abs(10^-7-yKr_h_bis2_vGB{1}));
    [a4 b4] = min(abs(10^-7-yKr_h_bis2_vGB{$}));
      

    subplot(4,4,xnf(ii))
    
        for i = 1:size(nn_BC,1)
            plot(xh_star_Kr_bis2_BC{nn_BC(i)},yKr_h_bis2_BC{nn_BC(i)},"color",id_color_BC(nn_BC(i)));
        end
    
        for i = 1:ntot_BC
            plot(xh_star_Kr_bis2_BC{i},yKr_h_bis2_BC{i},"color",id_color_BC(i));
        end
    
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$K_{r,\,BC}$","fontsize",xf2);
      
        xarrows([xh_star_Kr_bis2_BC{1}(b1);xh_star_Kr_bis2_BC{$}(b2)],[10^-7;10^-7],10*xar);
        
        xstring(8*10^5,1.5*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        ii = ii+1;
        
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [10^-3,10^-8;10^6,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lln";
        a.font_size = xf1;
//        atl = a.x_ticks.labels;
//        a.x_ticks.labels = "-"+atl;
        a.margins = xx_margins;
        a.auto_margins = "off";

        
    subplot(4,4,xnf(ii))

        for i = 1:ntot_vGB
            plot(xh_star_Kr_bis2_vGB{i},yKr_h_bis2_vGB{i},"color",id_color_vGB(i));
        end
    
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$K_{r,\,vGB}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
        
        a = gca();
        
        xarrows([xh_star_Kr_bis2_vGB{1}(b3);xh_star_Kr_bis2_vGB{$}(b4)],[10^-7;10^-7],10*xar);
        
        xstring(8*10^5,1.5*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        ii = ii+1;
        
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [10^-3,10^-8;10^6,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lln";
        a.font_size = xf1;
//        atl = a.x_ticks.labels;
//        a.x_ticks.labels = "-"+atl;
        a.margins = xx_margins;
        a.auto_margins = "off";


// position of the arrows for K(h)
    
    [a1 b1] = min(abs(10^-7-yKr_h_bis2_vGM{5}));
    [a2 b2] = min(abs(10^-7-yKr_h_bis2_vGM{$}));
    [a3 b3] = min(abs(10^-10-yKr_h_bis2_KG{12}));
    [a4 b4] = min(abs(10^-10-yKr_h_bis2_KG{3}));        // nn_KG_min:nn_KG_max

    subplot(4,4,xnf(ii))
        
        for i = 1:ntot_vGM
            plot(xh_star_Kr_bis2_vGM{i},yKr_h_bis2_vGM{i},"color",id_color_vGM(i));
        end
    
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$K_{r,\,vGM}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
        
        xarrows([xh_star_Kr_bis2_vGM{5}(b1);xh_star_Kr_bis2_vGM{$}(b2)],[10^-7;10^-7],10*xar);
        xarrows([10^-1;10^-1],[Kr_vGM(Se_hstar_vGM(10^-1,n_vGM(1)),n_vGM(1),l_vGM(1)); ..
        Kr_vGM(Se_hstar_vGM(10^-1,n_vGM($)),n_vGM($),l_vGM($))],10*xar);
        
        xstring(8*10^5,1.5*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        ii = ii+1;
        
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [10^-3,10^-8;10^6,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lln";
        a.font_size = xf1;
//        atl = a.x_ticks.labels;
//        a.x_ticks.labels = "-"+atl;
        a.margins = xx_margins;
        a.auto_margins = "off";

       
    subplot(4,4,xnf(ii))
    
        for i = nn_KG_min:nn_KG_max
            plot(xh_star_Kr_bis2_KG{nn_KG_max-i+nn_KG_min},yKr_h_bis2_KG{nn_KG_max-i+nn_KG_min},"color",id_color_KG(nn_KG_max-i+nn_KG_min));
        end
    
        xlabel("$|h^*|$","fontsize",xf2);
        ylabel("$K_{r,\,KG}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
        
        xarrows([xh_star_Kr_bis2_KG{12}(b3);xh_star_Kr_bis2_KG{3}(b4)],[10^-10;10^-10],15*xar);
        xarrows([10^-1;10^-1],[Kr_KG(Se_hstar_KG(10^-1,sigma_KG(nn_KG_min)),sigma_KG(nn_KG_min),l_KG(nn_KG_min)); ..
        Kr_KG(Se_hstar_KG(10^-1,sigma_KG(nn_KG_max)),sigma_KG(nn_KG_max),l_KG(nn_KG_max))],15*xar);
        
        a = gca();
        xstring(8*10^5,0.5*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        ii = ii+1;
        
        a = gca();
        a.axes_reverse = ["on","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [10^-3,10^-12;10^6,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "lln";
        a.font_size = xf1;
//        atl = a.x_ticks.labels;
//        a.x_ticks.labels = "-"+atl;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
 /// plots of D(Se) ///////////////////////////////////////////// 
 
// position of the arrows for D(Se)


    for i = 1:size(lambda_BC,1)
           DD(i) = abs(D_etoile_BC(0.2,lambda_BC(i),eta_BC(i)));
    end
    [c1 d1] = max(DD);
    
    for i = 1:size(n_vGB,1)
           DD(i) = abs(D_etoile_vGB(0.2,n_vGB(i),eta_vGB(i)));
    end
    [c2 d2] = max(DD);

    subplot(4,4,xnf(ii))
    
        for i = 1:ntot_BC
            plot(xSe_Dstar_bis2_BC{i}(1:($-1)),yDstar_Se_bis2_BC{i}(1:($-1)),"color",id_color_BC(i));
        end
          
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$D^*_{BC}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
    
        xarrows([0.2;0.2],abs([D_etoile_BC(0.2,lambda_BC(1),eta_BC(1));D_etoile_BC(0.2,lambda_BC(d1),eta_BC(d1))]),10*xar);
        xarrows([0.25;0.25],abs([D_etoile_BC(0.25,lambda_BC(d1),eta_BC(d1));D_etoile_BC(0.25,lambda_BC($),eta_BC($))]),10*xar);
    
        xstring(0,1.5*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        ii = ii+1;
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [0,10^-8;1,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "nln";
        a.font_size = xf1;
//      
    subplot(4,4,xnf(ii))
    
        for i = 1:ntot_vGB
            plot(xSe_Dstar_bis2_vGB{i}(1:($-1)),yDstar_Se_bis2_vGB{i}(1:($-1)),"color",id_color_vGB(i));
        end
          
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$D^*_{vGB}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
    
        xarrows([0.2;0.2],abs([D_etoile_vGB(0.2,n_vGB(1),eta_vGB(1));D_etoile_vGB(0.2,n_vGB(d2),eta_vGB(d2))]),10*xar);
        xarrows([0.25;0.25],abs([D_etoile_vGB(0.25,n_vGB(d2),eta_vGB(d2));D_etoile_vGB(0.25,n_vGB($),eta_vGB($))]),10*xar);
    
        xstring(0,1.5*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        ii = ii+1;
    
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [0,10^-8;1,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "nln";
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";
        
  // positions of arrows
        
    for i = 1:size(n_vGM,1)
           DD(i) = abs(D_etoile_vGM(0.2,n_vGM(i),l_vGM(i)));
    end
    [c1 d1] = max(DD);
    
    for i = nn_KG_min:nn_KG_max
           DD(i) = abs(D_etoile_KG(0.2,sigma_KG(i),l_KG(i)));
    end
    [c2 d2] = max(DD);
    
    subplot(4,4,xnf(ii))
    
        for i = 1:ntot_vGM
            plot(xSe_Dstar_bis2_vGM{i}(1:($-1)),yDstar_Se_bis2_vGM{i}(1:($-1)),"color",id_color_vGM(i));
        end
          
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$D^*_{vGM}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
    
        xarrows([0.2;0.2],abs([D_etoile_vGM(0.2,n_vGM(3),l_vGM(3));D_etoile_vGM(0.2,n_vGM(d1),l_vGM(d1))]),10*xar);
        xarrows([0.25;0.25],abs([D_etoile_vGM(0.25,n_vGM(d1),l_vGM(d1));D_etoile_vGM(0.25,n_vGM($),l_vGM($))]),7.5*xar);
    
        a = gca();
    
        xstring(0,1.5*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        ii = ii+1;
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [0,10^-8;1,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "nln";
        a.font_size = xf1;
        

    subplot(4,4,xnf(ii))
    
    n_min = 3;
    n_max = 17;
    
        for i = n_min:n_max
            plot(xSe_Dstar_bis2_KG{n_max-i+n_min}(1:($-1)),yDstar_Se_bis2_KG{n_max-i+n_min}(1:($-1)),"color",id_color_KG(n_max-i+n_min));
        end
          
        xlabel("$S_e$","fontsize",xf2);
        ylabel("$D^*_{KG}$","fontsize",xf2);
    //      title("Hydraulic conductivity curve - algo_2");
    
        xarrows([0.2;0.2],abs([D_etoile_KG(0.2,sigma_KG(nn_KG_min),l_KG(nn_KG_min));D_etoile_KG(0.2,sigma_KG(d2),l_KG(d2))]),12*xar);
        xarrows([0.25;0.25],abs([D_etoile_KG(0.25,sigma_KG(d2),l_KG(d2));D_etoile_KG(0.25,sigma_KG(nn_KG_max),l_KG(nn_KG_max))]),12*xar);
    
        a = gca();
        
        xstring(0,0.5*ry_abt2,alphabet(xnf(ii)));
        t = gce();                               // get the handle of the newly created object
        t.font_size=xfs_alphabet;
        ii = ii+1;
        
        a = gca();
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.data_bounds = [0,10^-12;1,10];
        a.tight_limits = ["on","on","off"];
        a.log_flags = "nln";
        a.font_size = xf1;
        a.margins = xx_margins;
        a.auto_margins = "off";

////////////////////// Savings images //////////////////////////////
//    
//    names = dir()(2);                       // 
//    
//    // savings en pdf ///
//    
//    name_pict_pdf(nfig_WRHCFs) = "WRHCFs.pdf";
//    name_pict_pdf(nfig_WRHCFs_xr) = "WRHCFs_xr.pdf";
//    
//    for i = 1:size(name_pict_pdf,1)
//        ltxt_pdf(i) = max(length(strstr(names,name_pict_pdf(i))));
//    end
//    
//    for i = 1:size(name_pict_pdf,1)
//        if ltxt_pdf(i) > 0 then,          // to refresh files
//            deletefile(name_pict_pdf(i));
//        end
//    end
//    
//    xs2pdf(nfig_WRHCFs,name_pict_pdf(nfig_WRHCFs));
//    xs2pdf(nfig_WRHCFs_xr,name_pict_pdf(nfig_WRHCFs_xr));
//    
//    
//    // savings en svg ////
//    
//    name_pict_svg(nfig_WRHCFs) = "WRHCFs.svg";
//    name_pict_svg(nfig_WRHCFs_xr) = "WRHCFs_xr.svg";
//    
//    for i = 1:size(name_pict_svg,1)
//        ltxt_svg(i) = max(length(strstr(names,name_pict_svg(i))));
//    end
//    
//    for i = 1:size(name_pict_svg,1)
//        if ltxt_svg(i) > 0 then,          // to refresh files
//            deletefile(name_pict_svg(i));
//        end
//    end
//    
//    xs2svg(nfig_WRHCFs,name_pict_svg(nfig_WRHCFs));
//    xs2svg(nfig_WRHCFs_xr,name_pict_svg(nfig_WRHCFs_xr));
    
//    
    disp("Computation time: " + string(floor(toc()/60*100)/100)+ " min");
