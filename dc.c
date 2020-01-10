#include <stdio.h>
#include <conio.h>
#include <math.h>
#define pi 3.141592654
#define muo (0.000001257)


int main(void)
{
    
    
                                                                 /*MAIN DIMENSION*/
    

    float P,n,eta,Bav,ac,Ki,p,Pa,Co,product,shi,val1,L,y,Li,f,t,tow,b,Va,D,Lmax,Dmax,Ddesired;
    double Dd;
    int nd;
    FILE *fp;
    fp=fopen("DCm.txt","w");
    fprintf(fp, "\n\n\n");
    
    printf("-------------------\n");
    printf("MAIN DIMENSION\n");
    printf("-------------------\n");
    fprintf(fp,"-------------------\n");
    fprintf(fp,"MAIN DIMENSION\n");
    fprintf(fp,"-------------------\n");

    printf("\nEnter the Power of the Motor(KW): ");
    scanf("%f",&P);
    printf("\nEnter the synchronous speed of the motor (RPM): ");
    scanf("%f",&n);
    printf("\nEnter the efficiency of the motor: ");
    scanf("%f",&eta);
    printf("\nEnter the value of magnetic loading (WB/m^2): ");
    scanf("%f",&Bav);
    printf("\nEnter the value of electric loading (A/m): ");
    scanf("%f",&ac);
    printf("\nEnter the value of stacking factor: ");
    scanf("%f",&Ki);
    printf("\nEnter the number of poles: ");
    scanf("%f",&p);
    if (P>50)
    {
        Pa=P;
    }
    if (P<50)
    {
        Pa=((1+2*eta)/(3*eta))*P;
    }
    printf("Pa=%f",Pa);
    Co=(pi*pi*Bav*ac)/1000;
    product=Pa/(Co*n/60);
    printf("\nD^2*L=%f",product);
    printf("\nEnter the value of pole arc to pole face ratio (.64 to .72): ");
    scanf("%f",&shi);
    // printf("enter the desired diameter: ");
    // scanf("%f",&Ddesired);

    val1=(product*p)/(shi*pi);

    Dd=cbrt(val1);
    printf("\nD = %f, val1 = %f",Dd,val1);   
    printf("\nThe calculated diameter is %f m",Dd);
    printf("\nEnter the standradized diameter: ");
    scanf("%f",&D);
    printf("\nD=%f",D);
    L=product/(D*D);
    printf("\nThe calculated length is %f m",L);
    printf("\nEnter the standradized length: ");
    scanf("%f",&L);
    tow=(pi*D)/p;
    b=shi*tow;
    Va=(pi*D*n)/60;
    printf("\nConsidering the length of onr ventilating duct is .13m");
    printf("\nEnter the width of ventilating duct: ");
    scanf("%f",&y);
    nd= L/.13;
    Li=Ki*(L-nd*y);
    f=(p*n)/(2*60);
    t=0.35;
    Lmax=7.5/(Bav*Va);
    Dmax=p/(pi*ac*Bav*L*Va)*1000;
    printf("the Lmax=%f         Dmax=%f",Lmax,Dmax);

                                                                // DESIGN SHEET

    fprintf(fp,"The power of the motor %f (KW)",P);
    fprintf(fp,"\nThe synchronous speed %f RPM is",n);
    fprintf(fp,"\nThe efficiency of the motor is %f",eta);
    fprintf(fp,"\nThe value of magnetic loading is %f WB/m^2",Bav);
    fprintf(fp,"\nThe value of electric loading is %f A/m",ac);
    fprintf(fp,"\nThe value of stacking factor is %f",Ki);
    fprintf(fp,"\nThe number of poles",p);
    fprintf(fp,"\nPower developed by armature is %f kW",Pa);
    fprintf(fp,"\nThe calculated diameter of armature is %f m",D);
    fprintf(fp,"\nThe D^2*L is %f m^3",product);
    fprintf(fp,"\nRatio of pole arc to pole face is %f m",D);
    fprintf(fp,"\nThe calculated diameter of armature is %f m",D);
    fprintf(fp,"\nThe standard Diameter is %f m",D);
    fprintf(fp,"\nThe calculated length of armature is %f m",L);
    fprintf(fp,"\nThe standard length is %f m",L);
    fprintf(fp,"\nThe Lmax=%f ________________ Dmax=%f",Lmax,Dmax);
    fprintf(fp,"\nPole pitch is %f m",tow);
    fprintf(fp,"\nPole arc is %f m",b);
    fprintf(fp,"\nPeripheral speed is %f m/s",Va);
    fprintf(fp,"\nConsidering the lingth of onr ventilating duct is .13m");
    fprintf(fp,"\nThe width of the ventilating duct is %f m",y);            
    fprintf(fp,"\nThe number of ventilating duct is %f",nd);
    fprintf(fp,"\nThe net iron length %f m",Li);
    fprintf(fp,"\nThe frequency of  %f m",f);
    fprintf(fp,"\nThe thickness of lamination used for machine is %f mm",t);



                                                         /*ARMATURE WINDING*/

    printf("\n-------------------\n");
    printf("ARMATURE WINDING\n");
    printf("-------------------\n");
    fprintf(fp,"-------------------\n");
    fprintf(fp,"ARMATURE WINDING\n");
    fprintf(fp,"-------------------\n");

    float V,II,val2,y1,If,y2,Ei,Ia,E,y3,la1,phi,Z,pk;
    printf("\nEnter the terminal voltage(V): ");
    scanf("%f",&V);
    II=(P*1000)/(eta*V);
    val2=Pa*n;
    printf("\nEnter the ratio of Line current to Field current-(reff,--fig.9.35): ");
    scanf("%f",&y1);
    If=(y1*II)/100;
    printf("\nEnter ratio of terminal voltage to internal voltage drop (as percentage)-(reff.--fig.9.35): ");
    scanf("%f",&y2);
    Ei=(y2*V)/100;
    Ia=II-If;
    E=V-Ei;
    printf("\nWe consider here the simplex wave winding: ");
    printf("\nEnter the type of winding   for WAVE TYPE '1' , for LAP TYPE '2' : ");
    scanf("%f",&pk);
    if (pk=1)
        {
            y3=p;
            fprintf(fp,"\nWe consider here the lap winding: ");
        }
    else if (pk=2)
    {
        printf("\nEnter the number of parallel path: ");
        scanf("%f",&y3);
        fprintf(fp,"\nWe consider here the simplex wave winding: ");
    }

    la1=Ia/y3;
    phi=Bav*tow*L;
    if((((E*y3)/(phi*(n/60)*p))-(((E*y3)/(phi*(n/60)*p))))>=0.5)
         Z=1+((E*y3)/(phi*(n/60)*p));
    else
         Z=((E*y3)/(phi*(n/60)*p));
    
                                                                 // DESIGN SHEET


    fprintf(fp,"\nThe terminal voltage is %f V:",V);
    fprintf(fp,"\nLine current is %f A",II);
    fprintf(fp,"\nProduct Pa*rpm is %f",val2);
    fprintf(fp,"\nThe ratio of line current to field current is %f percent",y1);
    fprintf(fp,"\nField current is %f A",If);
    fprintf(fp,"\nRatio of terminal voltage to internal voltage drop is %f",y2);
    fprintf(fp,"\nInternal voltage drop i %f Volt",Ei);
    fprintf(fp,"\nArmature current is %f A",Ia);
    fprintf(fp,"\nGenerated emf is %f Volt",E);
    fprintf(fp,"\nThe number of parallel path i %f",y3);
    fprintf(fp,"\nCurrent per parallel path i %f A",la1);
    fprintf(fp,"\nFlux per pole i %f Wb",phi);
    fprintf(fp,"\nNumber of armature conductor i %f",Z);



//                                                                 /*NUMBER OF SLOTS*/
    
    printf("-------------------\n");
    printf("NUMBER OF SLOTS\n");
    printf("-------------------\n");
    fprintf(fp,"-------------------\n");
    fprintf(fp,"NUMBER OF SLOTS\n");
    fprintf(fp,"-------------------\n");
  
    float SI,Su,S,Tc,u,C,sp;
    int Zc,Zs;
    noS:
    SI = (pi*D*1000)/35;
    Su=(pi*D*1000)/25;
    printf("\nNumber of slot lies between %f to %f",SI,Su);
    printf("\nEnter the Number of slots: ");
    scanf("%f",&S);
    printf("\nEnter the No of Turn per Coil: ");
    scanf("%f",&Tc);
    printf("\nEnter the Coil side per slot: ");
    scanf("%f",&u);
    C=(u*S)/2;
    Zc=(2*Tc*C);
    Zs=(Zc/S);
    sp=(pi*D)/S;
    phi=(phi*Z)/Zc;
    Bav=phi/(tow*L);
    printf("\nBav=%f\n",Bav);
    ac=(la1*Zc)/(pi*D);
    
                                                                    // DESIGN SHEET
    
    fprintf(fp,"\number of slot beween is %f and %f",SI,Su);
    fprintf(fp,"\number of slot is %f",S);
    fprintf(fp,"\number of turn per coil is %f",Tc);
    fprintf(fp,"\nColi size per slot is %f",u);
    fprintf(fp,"\nTotal number of cells is %f ",C);
    fprintf(fp,"\nTotal Number of armature conductor is %f",Zc);
    fprintf(fp,"\number of conductor per slot is %f", Zs);
    fprintf(fp,"\nin Slot pitch is %f m",sp);
    fprintf(fp,"\nModified value of flux per pole is %f WB",phi);
    fprintf(fp,"\nModified value of Magnetic loading is %f Wb/m^2", Bav);
    fprintf(fp,"\nModified value of Electric Loading is %f A/m",ac);




                                                                     /*CHECKS*/

    printf("-------------------\n");
    printf("CHECKS\n");
    printf("-------------------\n");
    fprintf(fp,"-------------------\n");
    fprintf(fp,"CHECKS\n");
    fprintf(fp,"-------------------\n");
 
    float Dc;
    int beta_c;
    if ((la1)<= 1500)
{
    fprintf(fp,"\nAmpere conductor per slot is lower than 1500");
    fprintf(fp,"\nValue of ampere conductor is %f A",la1*Zs);
    Dc=0.65*D;
    fprintf(fp,"\nDiameter of commutator is %f m",Dc);
    beta_c=(pi*Dc*10001/C);
    fprintf(fp,"\nValue of commutator segment is %fm",beta_c);
    if (beta_c<4)
    {
            printf("\nPitch of commutator segment is greater than 4.");
            fprintf(fp,"\nPitch of commutator segment is greater than 4.");
            printf("\nDo the require correction.\nStart new calculation for NO.of Slots");
            goto noS;
    }
}    
    else
    {
            printf("\nin Ampere conductor per slot is greater than 1500 A.");
            fprintf(fp,"\nAmpere conductor per slot is greater than 1500 A. It is %f A",la1*Zs);
            printf("\nDo the require correction.\nStart new calculation for NO.of Slots");
            goto noS;                                    
}

                                                                        // WINDING LAYOUT

    printf("-------------------\n");
    printf("WINDING LAYOUT\n");
    printf("-------------------\n");
    fprintf(fp,"-------------------\n");
    fprintf(fp,"WINDING LAYOUT\n");
    fprintf(fp,"-------------------\n");
 

    float yb,yc,Y,Yf;
    printf("\nEnter the Back pitch((2*C)/p+/-K) such that (yb-1)/u=integer: ");
    printf("\nWhere C=%0.0f;p=%0.0f;u=%0.Of;\n",C,p,u);
    printf("\nEnter the value of yb: ");
    scanf("%f",&yb);
    printf("\nEnter the Commutator pitch (2*(C+/-1)1/p) and such that (yc-1)/u=integer");
    printf("\nWhere C=%0.0f;p=%0.0f;u=%0.0f;\n",C,p,u);
    printf("\nEnter the value of yc: ");
    scanf("%f",&yc);
    Y=y3*yc;
    Yf=Y-yb;

                                                                         // DESIGN SHEET

    fprintf(fp,"\nThe Back pitch((2*C)/p+/-K) such that (yb-1)/u=integer");
    fprintf(fp, "\nBack pitch is %f",yb);
    fprintf(fp,"\nEnter the Commutator pitch(2*(C +/-1)/p) and such that (yc-1)/u=integer");
    fprintf(fp,"\nCommutator pitch is %f",yc);
    fprintf(fp,"\n Total winding pitch is %f",Y);
    fprintf(fp,"\nFront pitch is %f",Yf);


                                                                        // DESIGN OF SLOTS

    printf("-------------------\n");
    printf("DESIGN OF SLOTS\n");
    printf("-------------------\n");
    fprintf(fp,"-------------------\n");
    fprintf(fp,"DESIGN OF SLOTS\n");
    fprintf(fp,"-------------------\n");

    float delta_a,Aac,dm1,dm2,i,dmi1,dmi2,Ws,ds;
    printf("\nEnter the current Density in Armature winding A/mm^2 (Large strap wound armatures with very good normal ventilation: 4.5): ");    
    scanf("%f",&delta_a);
    Aac=la1/delta_a;
    printf("\nThe area of the armature conductor is %f (mm^2)",Aac);
    printf("\nReferring to table 17.1, select a copper conductor: ");
    scanf("%f%f",&dm1,&dm2);
    printf("\nEnter the insulation of the conductor(mm): ");
    scanf("%f",&i);
    dmi1=dm1+i;
    dmi2=dm2+i;
    printf("\nEnter the Total slot width (in m): ");
    scanf("%f",&Ws);
    printf("\nEnter the total slot depth (in m): ");
    scanf("%f",&ds);

                                                                         // DESIGN SHEET

    fprintf(fp,"\nThe current density in the armature core is %f A/mm^2",delta_a);
    fprintf(fp,"\nArea of the armature conductor is %f mm^2",Aac);
    fprintf(fp,"\nThe dimension of copper conductor is %f*%f",dm1,dm2);
    fprintf(fp,"\nInsulation of the conductor is %f mm",i);
    fprintf(fp,"\nInsulated dimension of the conductor is %f*%f",dmi1 ,dmi2);
    fprintf(fp,"\nTotal slot width is %f m",Ws);
    fprintf(fp,"\nTotal slot depth is %f m",ds);


                                                                        // CHECK THE FLUX DENSITY

    printf("-------------------\n");
    printf("CHECK THE FLUX DENSITY\n");
    printf("-------------------\n");
    fprintf(fp,"-------------------\n");
    fprintf(fp,"CHECK THE FLUX DENSITY\n");
    fprintf(fp,"-------------------\n");

    float ys1_3,Wt1_3,Bt1_3;
    ys1_3=pi*((D*1000)-((4*ds)/3))/S;
    Wt1_3=ys1_3-Ws*1000;
    Bt1_3=(p*phi*1000)/(shi*S*Li*Wt1_3);

                                                                        // DESIGN SHEET
    fprintf(fp,"\nSlot pitch at 1/3 height of the root is %f mm",ys1_3);
    fprintf(fp,"\nWidth of the tooth at 1/3 height from the root is %f mm",Wt1_3);
    fprintf(fp,"\nFlux density at 1/3 height of root is %f Wb/m^2",Bt1_3);


    //                                                                     // LENGTH OF AIR GAP

    // printf("-------------------\n");
    // printf("LENGTH OF AIR GAP\n");
    // printf("-------------------\n");
    // fprintf(fp,"-------------------\n");
    // fprintf(fp,"LENGTH OF AIR GAP\n");
    // fprintf(fp,"-------------------\n");


    // float ATa,ATg,Bg,Ig,kg;
    // ATa = (la1*Zc)/(2*p);
    // ATg=0.6*ATa;
    // Bg=Bav/shi;
    // printf("\nFlux Density in the air gap is %f Wb/m^2",Bg);
    // printf("\nEnter the Gap contraction factor: ");
    // scanf("%f",&kg);
    // Ig=ATg/(800000*Bg*kg);
    // printf("\nCalculated Length of air gap is %f m",Ig);
    // printf("\nEnter standard length of air gap in m: ");
    // scanf("%f",&Ig);

    //                                                                     // ARMATURE CORE

    // printf("-------------------\n");
    // printf("ARMATURE CORE\n");
    // printf("-------------------\n");
    // fprintf(fp,"-------------------\n");
    // fprintf(fp,"ARMATURE CORE\n");
    // fprintf(fp,"-------------------\n");

    // float phi_c,Bc,Ac,dc,Di;
    // phi_c=phi/2;
    // printf("\nEnter the density in the armature core (Wb/m^2): ");
    // scanf("%f",&Bc);
    // Ac=phi_c/Bc;
    // dc=Ac/Li;
    // printf("\nCalculated depth of armature core is %f m",dc);
    // printf("\nEnter the standard depth of armature core in m: ");
    // scanf("%f",&dc);
    // Ac=Li*dc;
    // Bc=phi_c/Ac;
    // Di=D-2*(dc+ds);

    //                                                                     // POLE SECTION

    // printf("-------------------\n");
    // printf("POLE SELECTION\n");
    // printf("-------------------\n");
    // fprintf(fp,"-------------------\n");
    // fprintf(fp,"POLE SELECTION\n");
    // fprintf(fp,"-------------------\n");


    // float Cl,phi_p,Bp,Ap,Lp,Lpi,bp;
    // printf("\nReferring table 9.7 Enter the value leakage coefficient of main poles: ");
    // scanf("%f",&Cl);
    // phi_p=Cl*phi;
    // printf("\nEnter the flux density in the pole body: ");
    // scanf("%f",&Bp);
    // Ap=phi_p/Bp;
    // Lp=L;
    // Lpi=Ki*L;
    // bp=Ap/Lpi;

    //                                                             // TENTATIVE DESIGN OF FIELD WINDING

    // printf("-------------------\n");
    // printf("TENTATIVE DESIGN OF FIELD WINDING\n");
    // printf("-------------------\n");
    // fprintf(fp,"-------------------\n");
    // fprintf(fp,"TENTATIVE DESIGN OF FIELD WINDING\n");
    // fprintf(fp,"-------------------\n");


    // float ATfl,y4,df,Sf,qf,AT,hf,hpl;
    // printf("\nEnter the ratio of mmf at full load to mmf of armature as %(0.9 to 1.25): ");
    // scanf("%f",&y4);
    // ATfl=(ATa,y4);
    // printf("\nReferring table 9.8, Enter the depth of the Field Winding in m(40 to 55): ");
    // scanf("%f",&df);
    // printf("\nEnter the value of Copper space factor(0.4 to 0.75): ");
    // scanf("%f",&qf);
    // printf("\nEnter the value of loss dissipation (700W/m^2): ");
    // scanf("%f",&Sf);
    // AT=10000*pow((qf*df*Sf),0.5);
    // hf=ATfl/AT;
    // hpl=(hf*1000+20+10)/1000;

    //                                                                         // YOKE

    // printf("-------------------\n");
    // printf("YOKE\n");
    // printf("-------------------\n");
    // fprintf(fp,"-------------------\n");
    // fprintf(fp,"YOKE\n");
    // fprintf(fp,"-------------------\n");

    // float phi_y,By,Ay,dy,Dy;
    // phi_y=phi/2;
    // printf("\nEnter the flux density in Yoke(Wb/m^2): ");
    // scanf("%f",&By);
    // Ay=phi_y/By;
    // dy=Ay/Lpi;
    // printf("\nCalculate depth of Yoke is %f m",dy);
    // printf("\nEnter the standard depth of Yoke in m: ");
    // scanf("&f",&dy);
    // Dy=D+2*(Ig/1000+hpl+dy);

    //                                                                     // MAGNETIC CIRCUIT


    // printf("\n-------------------\n");
    // printf("MAGNETIC CIRCUIT\n");
    // printf("-------------------\n");
    // fprintf(fp,"-------------------\n");
    // fprintf(fp,"MAGNETIC CIRCUIT\n");
    // fprintf(fp,"-------------------\n");


    // float val3,Kcs,Kgs,val4,Kcd,Kgd,Kg,Ks;
    // float att,ATt,Ic,atc,ATc,atp,ATp,ly,aty,ATy,ATf;
    // val3=Ws/(Ig);
    // printf("\nRatio of slot opening to gap length is %f",val3);
    // printf("\nReferring fig. 4.9 Enter the carter's coefficient corresponding to the ratio: ");
    // scanf("%f",&Kcs);
    // Kgs=(sp*1000)/((sp*1000)-(Kcs*(Ws*1000)));
    // val4=y/(Ig);
    // printf("\nRatio of duct width to gap length is %f",val4);
    // printf("\nReferring fig. 4.9 enter the carter's coefficient corresponding to the ratio: ");
    // scanf("%f",&Kcd);
    // Kgd=(L*1000)/((L*1000)-(nd*Kcd*(y*1000)));
    // Kg=Kgs*Kgd;
    // ATg=800000*Bg*Kg*(Ig/1000);
    // Ks=(L*1000*ys1_3)/(Li*1000*Wt1_3);
    // printf("\nThe value of flux density at 1/3 helght of Teeth is %f Wb/m^2",Bt1_3);
    // printf("\nEnter the value of mmf at teeth (A/m): ");
    // scanf("%f",&att);
    // ATt=att*ds;
    // Ic=(pi*(D-2*ds-dc))/(2*p);
    // printf("\nThe value of flux density at 1/3 helght of Core is Wb/m^2",Bc);
    // printf("\nEnter the value of mmf at Core(A/m): ");
    // scanf("%f",&atc);
    // ATc=atc*Ic;
    // printf("\nThe value of flux density at 1/3 height of Pole Body is %f Wb/m^2",Bp);
    // printf("\nEnter the value of mmf at Pole Body(A/m): ");
    // scanf("%f",&atp);
    // ATp=atp*hpl;
    // ly=(pi*(D+(2*Ig/1000)+2*hpl+dy))/(2*p);
    // printf("\nThe value of flux density at 1/3 height of Yoke is %f Wb/m^2",By);
    // printf("\nEnter the value of mmf at Yoke(A/m): ");
    // scanf("%f",&aty);
    // ATy=aty*ly;
    // ATf=ATg+ATt+ATc+ATp+ATy;
    // ATfl=ATf*1.15;

    //                                                             // DESIGN OF FIELD WINDING

    // printf("-------------------\n");
    // printf("DESIGN OF FIELD WINDING\n");
    // printf("-------------------\n");
    // fprintf(fp,"-------------------\n");
    // fprintf(fp,"DESIGN OF FIELD WINDINGS\n");
    // fprintf(fp,"-------------------\n");
    

    // float Ef,Lmt,row,af,d,di,Rf,Qf,Sc,cfd,theta;
    // int Tf;
    // printf("\nVoltage across the shunt field is %f Volt",V);
    // printf("\nand 20 percent of this voltage kept for reverse speed control.\n");
    // Ef=(V-0.2*V)/p;
    // Lmt=2*(Lp+bp+2*df);
    // printf("\nEnter the resistivity of the conductor for Copper-0.0211: ");
    // scanf("%f",&row);
    // af=(ATfl*row*Lmt)/Ef;
    // printf("\nReferring the table 17.7 enter the diameter of the conductor(mm): ");
    // scanf("%f",&d);
    // printf("\nReferring the table 17.7 enter the diameter with covering of the conductor(mm): ");
    // scanf("%f",&di);
    // Sf=0.75*pow((d/di),2);
    // Tf=(Sf*df*hf*1000000/af);
    // Rf=(Tf*row*Lmt)/af;
    // If=Ef/Rf;
    // ATfl=If*Tf;
    // Qf=(If*If)*Rf;
    // Sc=2*Lmt*(hf+df);
    // cfd=0.16/(1+(0.1*Va));
    // theta=(Qf*cfd)/Sc;

//                                                             // DESIGN OF COMUTATOR


//     printf("-------------------\n");
//     printf("DESIGN OF COMUTATOR\n");
//     printf("-------------------\n");
//     fprintf(fp,"-------------------\n");
//     fprintf(fp,"DESIGN OF COMUTATOR\n");
//     fprintf(fp,"-------------------\n");


//     float Vc,lb,y5,lbc,delta_b,ab,val5,tb,Wb,Ab,Lc,Lo;
//     Vc=pi*Dc*(n/60);
//     lb=(Ia*2)/p;
//     brush:
//     printf("\nCurrent Per brush arm is %f A",lb);
//     printf("\nThe current In each brush should not be more than 70A");
//     printf("\nEnter the number of brush per arm: ");
//     scanf("%f",&y5);
//     lbc=lb/y5;
//     if (lbc>70)
// {
//         printf("\nCurrent per brush arm is %f A",lbc);
//         printf("\nChange the Number of brush");
//         getchar();
//         goto brush;
// }
//     delta_b=0.1;
//     ab=lbc/delta_b;
//     val5=2.5;
//     tb=val5*beta_c;
//     Wb=ab/tb;
//     printf("\nArea of each brush is %f mm^2",ab);
//     printf("\nThickness of brush is %0.3f mm and width of brush is %0.3f mm",tb,Wb);
//     printf("\nEnter the dimension of brush (in mm) referring the table: ");
//     scanf("%f%f",&tb,&Wb);
//     ab=(tb/1000)*(Wb/1000);
//     Ab=y5*ab;
//     Lc=(y5*(Wb+5)+10+10)/1000;
//     Lo=Lc+0.02;




//                                                                         // LOSSES


//     printf("-------------------\n");
//     printf("LOSSES\n");
//     printf("-------------------\n");
//     fprintf(fp,"-------------------\n");
//     fprintf(fp,"LOSSES\n");
//     fprintf(fp,"-------------------\n");


//     float Pbc,Vbd,val6,val7,Pbf,Lt,Bsc,cl,theta_l;
//     printf("\nEnter the brush contact drop(1V per brush): ");
//     scanf("%f",&Vbd);
//     Pbc=2*Vbd*Ia;
//     val6=20000;
//     val7=0.15;
//     Pbf=val6*val7*y5*Ab*Vc;
//     Lt=Pbc+Pbf;
//     Bsc=pi*Dc*Lc+0.01;
//     cl=0.025/(1+0.1*Vc);
//     theta_l=(Lt*cl)/Bsc;

//                                                                     // DESIGN OF INTER POLES

//     printf("-------------------\n");
//     printf("DESIGN OF INTER POLES\n");
//     printf("-------------------\n");
//     fprintf(fp,"-------------------\n");
//     fprintf(fp,"DESIGN OF INTER POLES\n");
//     fprintf(fp,"-------------------\n");


//     float Wc,Igi,Wip,lamda_t,lamda_s,Le,Loo,bo,lamda_o;
//     float lamda,tow_c,Erav,Em,Bgim,Atgi,Kgi,ATi,delta_i,ai;
//     float y6,y7;
//     int Ti;
//     Wc=((u/2-y3/2)*beta_c+tb)*(D/Dc);
//     Igi=1.2*Ig;
//     printf("\nMax Width of Interpole sholld be greater than 1.5 times of slot pitch");
//     printf("\nEnter the value of Width of Interpolet (in mm) such that the Inter pole");
//     printf("\nface will cover commutation zone and also approximate to slot pitch");
//     printf("\nThe value of slot pitch is %0.34f m: ",sp);
//     scanf("%f",&Wip);
//     printf("\nEnter the specific slot permeance: ");
//     scanf("%f",lamda_s);
//     lamda_t=(muo*Wip)/(6*Igi);
//     Le=0.3*tow+0.125*(ds/1000);
//     Loo=pow(((tow*tow/4)+(Le*Le)),0.5);
//     printf("\nEnter the peripheri of the coil side(in m): ");
//     scanf("%f",&bo);
//     lamda_o=0.000001*(Loo/L)*(0.23*log10(Loo/bo)+0.07);
//     lamda=(lamda_s+lamda_t+lamda_o);
//     tow_c=((u/2-y3/2)*beta_c+tb)/Vc;
//     Erav=4000*Tc*lamda*L*la1*(Zs/tow_c);
//     printf("\nEnter the maximum reactance voltage(Volt): ");
//     scanf("%f",&Em);
//     Bgim=Em/(L*Va);
//     printf("\nThe length of air gap under the interpole is greater than the main pole and");
//     printf("\ntherefore, the gap contraction factor for \nthe interpoles is smaller than Kg(%f))",Kg);
//     printf("\nEnter the value gap contraction factor for interpoles: ");
//     scanf("%f",&Kgi);
//     Atgi=800000*Bgim*Kgi*Igi/1000;
//     ATi=Atgi+ATa;
//     Ti=(ATi/Ia);
//     printf("\nEnter the current density for the interpoles(aA/mm^2): ");
//     scanf("%f",&delta_i);
//     ai=Ia*delta_i;
//     printf("\nArea of the inpterpole winding is %f mm^2",ai);
//     printf("\nEnter the standard dimension from the tablet(in mm): ");
//     scanf("%f%f",&y6,y7);

//                                                                     // LOSSES AND EFFICIENCY


//     printf("-------------------\n");
//     printf("LOSSES AND EFFICIENCY\n");
//     printf("-------------------\n");
//     fprintf(fp,"-------------------\n");
//     fprintf(fp,"LOSSES AND EFFICIENCY\n");
//     fprintf(fp,"-------------------\n");


//     float y8,Pbw,Pt,Wmt,Wt,Lspt,TI,Lspc,Closs,Tir;
//     float IRI,Lmta,ra,Lca,Lcs,Lmti,ri,Lci,Lbc,Loss,eff;
//     printf("\nReffering table 9.11 enter the bearing and\n windage losses (in percentage): ");
//     scanf("%f",&y8);
//     Pbw=(P*1000*y8)/100;
//     Pt=Pbf+Pbw;
//     Wmt=pi*(D-ds)/S-Ws;
//     Wt=S*L*Wmt*ds*7800;
//     Lspt=0.06*f*Bt1_3*Bt1_3+0.008*f*f*Bt1_3*Bt1_3*t*t;
//     TI=Wt*Lspt;
//     printf("\nEnter the weight of armature core(kg): ");
//     scanf("%f",&Wc);
//     Lspc=0.06*f*Bc*Bc+0.005*f*f*Bc*Bc*t*t;
//     Closs=Lspc*Wc;
//     Tir=Ti+Closs;
//     IRI=1.2*Tir;
//     Lmta=2*L+2.3*tow+5*ds;
//     printf("\nEnter the armature resistance (ohm): ");
//     scanf("%f",&ra);
//     Lca=Ia*Ia*ra;
//     Lcs=V*If;
//     printf("\nEnter the length of mean turn for Interpolation (in mm): ");
//     scanf("%f",&Lmti);
//     ri=y5*Ti*row*Lmti/ai;
//     Lci=Ia*Ia*ri;
//     Lbc=Pbc;
//     Loss=Lca+Lcs+Lci+Lbc+IRI+Pt;
//     eff=(P*1000*100)/((P*1000)+Loss);
//     printf("\n%f",eff);

//     printf("\n%f %f %f %f",Pt,Wmt,IRI,Lmta);

//     getchar();
//     return 0;


fclose(fp);
}