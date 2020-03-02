#include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #define PI 3.14159265359
    #define l1 1
    #define l2 1

    double I1=0,I2=0, g=9.81,theta1i=PI/2,theta2i=PI/2,p1i=0,p2i=0;
    static double System(int i,double theta1, double theta2, double p1, double p2, double m[6]);
    static double massdensityf(double i,double r);
    int main()
    {
    FILE *f1;FILE *f2;FILE *f3;
    f1 = fopen("variables.txt", "w");
    f2 = fopen("parameters.txt", "w");
    f3 = fopen("poincare.txt", "w");
    int i,j,c=0;
    double N=1000,Nm=10000,r1=0,r2=0, a=0,B=0,m1=0,m2=0,ma[6], angle,Ei;
    double h=1/N,hm1=l1/Nm,hm2=l2/Nm;

    //Center of mass and total mass calculation
    for (i=1;i<Nm;i++)
    {
        a+=0.5*massdensityf(1,r1+hm1)*(r1+hm1)*hm1+0.5*massdensityf(1,r1)*r1*hm1;
        B+=0.5*massdensityf(2,r2+hm2)*(r2+hm2)*hm2+0.5*massdensityf(2,r2)*r2*hm2;
        m1+=0.5*massdensityf(1,r1+hm1)*hm1+0.5*massdensityf(1,r1)*hm1;
        m2+=0.5*massdensityf(2,r2+hm2)*hm2+0.5*massdensityf(2,r2)*hm2;
        r1+=hm1;
        r2+=hm2;
    }
    a=a/(l1*m1);
    B=B/(l2*m2);
    r1=-a;
    r2=-B;
    for (i=1;i<Nm;i++)
    {
        I1+=0.5*massdensityf(1,r1+a+hm1)*(r1+hm1)*(r1+hm1)*hm1+0.5*massdensityf(1,r1+a)*r1*r1*hm1;
        I2+=0.5*massdensityf(2,r2+B+hm2)*(r2+hm2)*(r2+hm2)*hm2+0.5*massdensityf(2,r2+B)*r2*r2*hm2;
        r1+=hm1;
        r2+=hm2;
    }

    ma[0]=I1;ma[1]=I2;ma[2]=a;ma[3]=B;ma[4]=m1;ma[5]=m2;

    /*Runge Kutta*/

    double k[4][4],x[]={theta1i,theta2i,p1i,p2i},x1[]={theta1i,theta2i,p1i,p2i}, t=0,E, T=200,theta1dt,theta2dt,errorE=0;
    while(t<T)
    {
        theta1dt=((x[2]*(m2*l2*l2*B*B+I2)-x[3]*m2*l1*l2*B*cos(x[0]-x[1]))/((m1*l1*l1*a*a+m2*l1*l1+I1)*(m2*l2*l2*B*B+I2)-pow(m2*l1*l2*B*cos(x[0]-x[1]),2)));
        theta2dt=((x[3]*(m1*l1*l1*a*a+m2*l1*l1+I1)-x[2]*m2*l1*l2*B*cos(x[0]-x[1]))/((m2*l2*l2*B*B+I2)*(m1*l1*l1*a*a+m2*l1*l1+I1)-pow(m2*l1*l2*B*cos(x[0]-x[1]),2)));
        E=0.5*pow(theta1dt,2)*(m1*l1*l1*a*a+m2*l1*l1+I1)+0.5*pow(theta2dt,2)*(m2*l2*l2*B*B+I2)+m2*l1*l2*theta1dt*theta2dt*B*cos(x[0]-x[1])-m1*g*l1*a*cos(x[0])-m2*g*(l1*cos(x[0])+l2*B*cos(x[1]))+m1*g*l1*a+m2*g*(l1+l2*B);
        if (t==0){Ei=E;}
        if (errorE<fabs(E-Ei)/E){errorE=fabs(E-Ei)/E;}

        //printf("%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.5lf\t%8.3lf\t%e\n",t, x[0],x[1],x[2],x[3],theta1dt,theta2dt,E,Ei,fabs(E-Ei));
        if (c % 40==0 || c==0){fprintf(f1,"%8.31f\t%8.31f\t%8.3lf\t%8.3lf\t%8.3lf\t%8.5lf\t%8.3lf\t%e\n",t, x[0],x[1],x[2],x[3],E,Ei,fabs(E-Ei));}


        for (j=0;j<4;j++)
        {
            for (i=0;i<4;i++)
            {
                if (j==0){k[i][0]=h*System(i,x[0], x[1], x[2], x[3],ma);}
                if (j==1){k[i][1]=h*System(i,x[0]+k[0][0]/2.0,x[1]+k[1][0]/2.0,x[2]+k[2][0]/2.0,x[3]+k[3][0]/2.0,ma);}
                if (j==2){k[i][2]=h*System(i,x[0]+k[0][1]/2.0,x[1]+k[1][1]/2.0,x[2]+k[2][1]/2.0,x[3]+k[3][1]/2.0,ma);}
                if (j==3){k[i][3]=h*System(i,x[0]+k[0][2],x[1]+k[1][2],x[2]+k[2][2],x[3]+k[3][2],ma);}
            }
        }

        for (i=0;i<4;i++)
        {
            x1[i]=x[i];//Save x_(n-1)
            x[i]=x[i]+(k[i][0]+2.0*k[i][1]+2.0*k[i][2]+k[i][3])/6.0;
        }
        //Poincare section, when theta2=0
        if (cos(x[1])>0 && sin(x[1])*sin(x1[1])<0)
        {
              angle=x[0];
              if (angle>0){angle=angle-2*PI*(floor(angle/(2*PI)));}
              if (angle<0){angle=angle-2*PI*(ceil(angle/(2*PI)));}
              if (fabs(angle)>PI){angle=-(((angle>0)-(angle<0))*2*PI-angle);}

              printf("%8.3f\t%8.3f\t%8.3f\t%8.3f\n",t,angle);
              fprintf(f3,"%e\t%e\n",angle,theta1dt);
        }

        t=t+h;
        c+=1;

    }
    printf("\n\n%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%10.5e\t%10.5e\t%10.5e\n",m1, m2,I1,I2,a,B,E,Ei,errorE);
    fprintf(f2,"m1\t\tm2\t\tI1\t\tI2\t\ta\t\tB\t\tE\t\tEi\t\tRelative Energy Error\n");
    fprintf(f2,"%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%10.5e\t%10.5e\t%10.5e\n",m1, m2,I1,I2,a,B,E,Ei,errorE);
    fclose(f1);fclose(f2);fclose(f3);
    return(0);

}

double massdensityf(double i,double r)
{
    if (i==1){return 1/l1;}
    if (i==2){return 1/l2;}
}

double System(int i,double theta1, double theta2, double p1, double p2, double ma[])
    {
    double I1=ma[0],I2=ma[1],a=ma[2],B=ma[3],m1=ma[4],m2=ma[5];
    double theta1dt=((p1*(m2*l2*l2*B*B+I2)-p2*m2*l1*l2*B*cos(theta2-theta1))/((m1*l1*l1*a*a+m2*l1*l1+I1)*(m2*l2*l2*B*B+I2)-pow(m2*l1*l2*B*cos(theta2-theta1),2)));
    double theta2dt=((p2*(m1*l1*l1*a*a+m2*l1*l1+I1)-p1*m2*l1*l2*B*cos(theta2-theta1))/((m2*l2*l2*B*B+I2)*(m1*l1*l1*a*a+m2*l1*l1+I1)-pow(m2*l1*l2*B*cos(theta2-theta1),2)));
    if (i==0){return  (theta1dt);}
    if (i==1){return  (theta2dt);}
    if (i==2){return  (m2*l1*l2*B*theta1dt*theta2dt*sin(theta2-theta1)-g*sin(theta1)*l1*(m1*a+m2));}
    if (i==3){return  (-m2*l1*l2*B*theta1dt*theta2dt*sin(theta2-theta1)-l2*m2*g*B*sin(theta2));}
    }
