// description
    
    tic();  // clock
    
    // Computation of sorptivity and c_p parameter by L. Lassabatere et al. (2021)
    /// long computations --> please wait!!

// Initialization
  
    clear
    mode(0);
    ieee(1);
    clearglobal;
    i_fig = 1;
    
// Working repertory and loading functions
  
    T_wd = pwd();
    
    T_rep = T_wd + "\0_functions";
    
    chdir(T_rep);
    
    exec('0_general_tool_functions.sci',-1);
    exec('1_WRCHCF_func.sci',-1);
    
    close_windows(20);
    
    chdir(T_wd);
  
//// Entries: xh and Se, and shape index /////////////////////////////// 
    
    xSe0 = 0;
    xSef = 1;
    
    xR = linspace(0.01,0.99,99)'+10^-4;         
    // values of m or parameter between 0 and 1
    // 10^-4 was added to avoid the numerical indetermination for m = 1/2
    
//// Hydraulic shape parameters /////////////////////////////// 

 // BC
    
    xlambda_BC = 2*xR./(1-xR);
    eta_BC = 2./xlambda_BC+3;
    
 // vGB
    
    xn_vGB = floor(2./(1-xR)*100)/100;
    xm_vGB = 1-2./xn_vGB;
    eta_vGB = 2./(xn_vGB-2)+3;
    
 // vGM
    
    xn_vGM = floor(1./(1-xR)*100)/100;
    xm_vGM = 1-1./xn_vGM;
    l_vGM = 0.5*ones(xn_vGM);
    
 // KG
    
    xsigma_KG = 1./xR-1;
    l_KG = 0.5*ones(xsigma_KG);
    
/// Computation of dimensionless sorptivity^2 (c_p)      
    
 // BC
    
    for i = 1:size(xlambda_BC,1)                                       // equation S2_BC(psi) of the paper
        psi_bc(i) = xlambda_BC(i)/(xlambda_BC(i)+2);
        yS2_BC(i) = cp_function_BC_psi(psi_bc(i));
    end
    
 // vGB
    
    for i = 1:size(xm_vGB,1)
        yS2_vGB(i) = cp_function_vGB_m(xm_vGB(i));
    end
    
 // vGM
    
    for i = 1:size(xn_vGM,1)
        xm_vGM(i) = 1-1/xn_vGM(i);
        yS2_vGM(i) = cp_function_vGM_m(xm_vGM(i),l_vGM(i));
    end
    
 // KG
    
    for i = 17:size(xsigma_KG,1)
        yS2_KG(i) = abs(S_etoile_2_KG_comb(xSe0,xSef,xsigma_KG(i),l_KG(i)));
    end
    
/// Computation of dimensionless sorptivity (c_p^1/2)
    
    for i = 1:size(yS2_vGB,1)
        yS_vGB(i) = sqrt(yS2_vGB(i));
    end
    
    for i = 1:size(yS2_vGM,1)
        yS_vGM(i) = sqrt(yS2_vGM(i));
    end
    
    for i = 1:size(yS2_BC,1)
        yS_BC(i) = sqrt(yS2_BC(i));
    end
    
    for i = 1:size(yS2_KG,1)
        yS_KG(i) = sqrt(yS2_KG(i));
    end

    
///////////////////////////// Figure 3 of the paper /////////////////////////////////////////////////////////////

  // fig. options
  
    xf1 = 3;     // font_size
    xf2 = 4;     // font_size
    xleg = 4;
    xth = 1.;
    
    xx_bounds = [0,0;1,4];
    
  // subfigures legends
    
    x_abt = 0.02;
    y_abt = 3.55;
    
    alphabet = "$\bf {\large {("+["a" "b" "c" "d" "e" "f" "g" ..
                    "h" "i" "j" "k" "l" "m" "n" ..
                    "o" "p" "q" "r" "s" "t" "u" ..
                    "v" "w" "x" "y" "z"]+")}}$";
                    
    ii = 1;
     
  // fig. axes ticks
     
    x_t_location=[0:0.1:1];
    x_t_label='$'+string(x_t_location)+'$';
    ticks_x = tlist(["ticks","locations","labels"], x_t_location, x_t_label);
    
    y_t_location=[0:0.5:4];
    y_t_label='$'+string(y_t_location)+'$';
    ticks_y = tlist(["ticks","locations","labels"], y_t_location, y_t_label);
     

  // figs
    
    clf(i_fig);
    scf(i_fig);
    
    fig = gcf();
    fig.figure_size = [1050,900];
    
    i_fig = i_fig+1;
    
    subplot(221)
    
        plot([0;xR;1],[2;yS_BC;sqrt(2)].^2,"red -");
        plot(xR,2*ones(xR),"--");
        
        leg = legend("$c_p\left( x \right)$","$c_{p,d}=2$");
        leg.font_size = xleg;
        leg.line_mode = "off";
        leg.fill_mode = "off";
        
        xlabel("$x=\frac{\lambda_{BC}}{2+\lambda_{BC}}$","fontsize",xf2);
        ylabel("$c_{p,BC}\left( x \right)$","fontsize",xf2);
        
        a = gca();
        a.data_bounds = xx_bounds;
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nnn";
        a.font_size = xf1;
        a.auto_ticks=['off' 'off' 'on'];
        a.x_ticks=ticks_x;
        a.y_ticks=ticks_y;
        
        xstring(x_abt,y_abt,alphabet(ii));
        
        ii = ii + 1;
        
    subplot(222)
    
        plot([0;xR;1],[0;yS_KG;sqrt(2)].^2,"red -");
        plot(xR,2*ones(xR),"--");
        
        xlabel("$x=\frac{1}{1+\sigma_{KG}}$","fontsize",xf2);
        ylabel("$c_{p,KG}\left( x \right)$","fontsize",xf2);
        
        a = gca();
        a.data_bounds = xx_bounds;
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nnn";
        a.font_size = xf1;
        a.auto_ticks=['off' 'off' 'on'];
        a.x_ticks=ticks_x;
        a.y_ticks=ticks_y;
        
        xstring(x_abt,y_abt,alphabet(ii));
        
        ii = ii + 1;
        
        
    subplot(223)
    
        plot([0;xm_vGB;1],[gamma(1/2);yS_vGB;sqrt(2)].^2,"red -");
        plot(xm_vGB,2*ones(xm_vGB),"--");
        
        xlabel("$x=m_{vGB}$","fontsize",xf2);
        ylabel("$c_{p,vGB}\left( x \right)$","fontsize",xf2);
        
        a = gca();
        a.data_bounds = xx_bounds;
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nnn";
        a.font_size = xf1;
        a.auto_ticks=['off' 'off' 'on'];
        a.x_ticks=ticks_x;
        a.y_ticks=ticks_y;
        
        xstring(x_abt,y_abt,alphabet(ii));
        
        ii = ii + 1;
        
    subplot(224)
    
        plot([0;xm_vGM;1],[0;yS_vGM;sqrt(2)].^2,"red -");
        plot(xm_vGM,2*ones(xm_vGM),"--");
        
        xlabel("$x=m_{vGM}$","fontsize",xf2);
        ylabel("$c_{p,vGM}\left( x \right)$","fontsize",xf2);
        
        a = gca();
        a.data_bounds = xx_bounds;
        a.axes_reverse = ["off","off","off"];
        a.axes_visible = ["on","on","on"];
        a.log_flags = "nnn";
        a.font_size = xf1;
        a.auto_ticks=['off' 'off' 'on'];
        a.x_ticks=ticks_x;
        a.y_ticks=ticks_y;
    
        xstring(x_abt,y_abt,alphabet(ii));
        
        ii = ii + 1;
    
//////// Saving pictures //////////////////////////////////
    
//    name_pict = "parameter_cp.pdf";
//    name_pict_2 = "parameter_cp.svg";
//    
//    names = dir()(2);
//    ltxt = length(strstr(names,name_pict));
//    ltxt2 = length(strstr(names,name_pict_2));
//    
//    if max(ltxt) > 0 then,          // to refresh files
//        deletefile(name_pict);
//    end
//    
//    if max(ltxt) > 0 then,          // to refresh files
//        deletefile(name_pict_2);
//    end
//    
//    xs2pdf(1,name_pict);
//    xs2svg(1,name_pict_2);
    
//////// SAVINGS //////////////////////////////////

    nc = [2:2:98];

    xR_res = [0;xR(nc);1];
    
    xlambda_BC_res = [0;xlambda_BC(nc);max(xlambda_BC)];
    xn_vGB_res = [2;xn_vGB(nc);max(xn_vGB)];
    xn_vGM_res = [1;xn_vGM(nc);max(xn_vGM)];
    xsigma_KG_res = [max(xsigma_KG);xsigma_KG(nc);min(xsigma_KG)];
    
    yS_BC_res = [2;yS_BC(nc);sqrt(2)].^2;
    yS_vGB_res = [gamma(1/2);yS_vGB(nc);sqrt(2)].^2;
    yS_vGM_res = [0;yS_vGM(nc);sqrt(2)].^2;
    yS_KG_res = [0;yS_KG(nc);sqrt(2)].^2;
    
    M1 = [xR_res xlambda_BC_res xn_vGB_res xn_vGM_res xsigma_KG_res];
    M1 = [xR_res];
    M2 = [yS_BC_res yS_vGB_res yS_vGM_res yS_KG_res];
    
    Mtot = [M1 M2];
    // values of cp
    // for Table 1
    
    disp("Computation time: " + string(floor(toc()/60*100)/100)+ " min");
