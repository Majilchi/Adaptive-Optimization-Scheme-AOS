
function f=BenFunctions(L,FunIndex,dim)
switch FunIndex
    
    case 1   %Stepint
        f=25+sum(floor(L));
        
    case 2    %Step
        f=sum(floor(L+0.5).^2);
        
    case 3   %Sphere;
        f=(sum(L.^2));
        
    case 4  %SumSquares
        f=0;
        for i=1:dim
            f=f+i*L(i)^2;
        end
        
    case 5  %Quartic
        dim=length(L);
        f=sum((1:dim).*(L.^4))+rand;
        
    case 6    %Beale
        f1=(1.5-L(1)+L(1)*L(2))^2;
        f2=(2.25-L(1)+L(1)*L(2)^2)^2;
        f3=(2.625-L(1)+L(1)*L(2)^3)^2;
        f=f1+f2+f3;
        
    case 7   %Easom
        f1=-cos(L(1))*cos(L(2));
        f2=exp(-(L(1)-pi)^2-(L(2)-pi)^2);
        f=f1*f2;
        
    case 8  %Matyas
        f1=0.26*(L(1)^2+L(2)^2);
        f2=-0.48*L(1)*L(2);
        f=f1+f2;
        
    case 9  %Colville
        f1=100*(L(1)^2-L(2))^2;
        f2=(L(1)-1)^2;
        f3=(L(3)-1)^2;
        f4=90*(L(3)^2-L(4))^2;
        f5=10.1*((L(2)-1)^2+(L(4)-1)^2);
        f6=19.8*(L(2)-1)*(L(4)-1);
        f=f1+f2+f3+f4+f5+f6;
        
    case 10  %Trid6
        f1=(L(1)-1)^2;
        f2=0;
        for i=2:dim
            f1=f1+(L(i)-1)^2;
            f2=f2+L(i)*L(i-1);
        end
        f=f1-f2;
        
    case 11 %Trid10
        f1=(L(1)-1)^2;
        f2=0;
        for i=2:dim
            f1=f1+(L(i)-1)^2;
            f2=f2+L(i)*L(i-1);
        end
        f=f1-f2;
        
    case 12  %Zakharov
        f1=0;
        f2=0;
        for i=1:dim
            f1=f1+L(i)^2;
            f2=f2+0.5*i*L(i);
        end
        f=f1+f2^2+f2^4;
        
    case 13  %Powell
        f=0;
        for i=1:(dim/4)
            f1=(L(4*i-3)+10*L(4*i-2))^2;
            f2=5*(L(4*i-1)-L(4*i))^2;
            f3=(L(4*i-2)-2*L(4*i-1))^4;
            f4=10*(L(4*i-3)-L(4*i))^4;
            f=f+f1+f2+f3+f4;
        end
        
    case 14 %Schwefel 2.22
        dim=length(L);
        f1=0;
        f2=1;
        for i=1:dim
            f1 =f1+ abs(L(i));
            f2 = f2*abs(L(i));
        end
        f=f1+f2;
        
    case 15 %Schwefel 1.2
        dim=length(L);
        f=0;
        for i=1:dim
            f=f+sum(L(1:i))^2;
        end
        
    case 16 %Rosenbrock
        dim=length(L);
        f=sum(100*(L(2:dim)-(L(1:dim-1).^2)).^2+(L(1:dim-1)-1).^2);
        
    case 17  %Dixon-Price
        dim=length(L);
        f0=(L(1)-1)^2;
        f1=0;
        for i=2:dim
            y=i*(2*L(i)^2-L(i-1))^2;
            f1=f1+y;
        end
        f=f0+f1;
        
    case 18  %Foxholes
        aij=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;...
            -32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];
        for i=1:25
            s(i)=sum((L'-aij(:,i)).^6);
        end
        f=(1/500+sum(1./((1:25)+s))).^(-1);
        
    case 19 %Branin
        f=(L(2)-(L(1)^2)*5.1/(4*(pi^2))+5/pi*L(1)-6)^2+10*(1-1/(8*pi))*cos(L(1))+10;
        
    case 20 %Bohachevsky1
        f1=L(1)^2;
        f2=2*L(2)^2;
        f3=-0.3*cos(3*pi*L(1));
        f4=-0.4*cos(4*pi*L(2));
        f=f1+f2+f3+f4+0.7;
        
    case 21 %Booth
        f1=(L(1)+2*L(2)-7)^2;
        f2=(2*L(1)+L(2)-5)^2;
        f=f1+f2;
        
    case 22  %Rastrigin
        dim=length(L);
        f=sum(L.^2-10*cos(2*pi.*L))+10*dim;
        
    case 23  %Schwefel
        f=sum(-L.*sin(sqrt(abs(L))));
        
    case 24  %Michalewicz2
        dim=length(L);
        a=10;
        f1=0;
        for i=1:dim
            f1=f1+sin(L(i))*(sin(i*L(i)^2/pi))^(2*a);
        end
        f=-f1;
        
    case 25 %Michalewicz5
        dim=length(L);
        a=10;
        f1=0;
        for i=1:dim
            f1=f1+sin(L(i))*(sin(i*L(i)^2/pi))^(2*a);
        end
        f=-f1;
        
    case 26   %Michalewicz10
        dim=length(L);
        a=10;
        f1=0;
        for i=1:dim
            f1=f1+sin(L(i))*(sin(i*L(i)^2/pi))^(2*a);
        end
        f=-f1;
        
    case 27  %Schaffer
        f1=(sin(L(1)^2-L(2)^2))^2-0.5;
        f2=(1+0.001*(L(1)^2+L(2)^2))^2;
        f=0.5+f1/f2;
        
    case 28  %Six Hump Camel
        f=4*(L(1)^2)-2.1*(L(1)^4)+(L(1)^6)/3+L(1)*L(2)-4*(L(2)^2)+4*(L(2)^4);
        
    case 29  %Bohachevsky2
        f1=L(1)^2;
        f2=2*L(2)^2;
        f3=-0.3*cos(3*pi*L(1))*cos(4*pi*L(2));
        f=f1+f2+f3+0.3;
        
    case 30   %Bohachevsky3
        f1=L(1)^2;
        f2=2*L(2)^2;
        f3=-0.3*cos(3*pi*L(1)+4*pi*L(2));
        f=f1+f2+f3+0.3;
        
    case 31  %Shubert
        f1=0;
        f2=0;
        for i=1:5
            s1=i*cos((i+1)*L(1)+i);
            s2=i*cos((i+1)*L(2)+i);
            f1=f1+s1;
            f2=f2+s2;
        end
        f=f1*f2;
        
    case 32   %GoldStein-Price
        f=(1+(L(1)+L(2)+1)^2*(19-14*L(1)+3*(L(1)^2)-14*L(2)+6*L(1)*L(2)+3*L(2)^2))*...
            (30+(2*L(1)-3*L(2))^2*(18-32*L(1)+12*(L(1)^2)+48*L(2)-36*L(1)*L(2)+27*(L(2)^2)));
        
    case 33  %Kowalik
        a=[0.1957 0.1947 0.1735 0.16 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
        b=1./[0.25 0.5 1 2 4 6 8 10 12 14 16];
        f=sum((a-((L(1).*(b.^2+L(2).*b))./(b.^2+L(3).*b+L(4)))).^2);
        
    case 34  %Shekel5
        a=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
        c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
        f=0;
        for i=1:5
            f=f-1/((L-a(i,:))*(L-a(i,:))'+c(i));
        end
        
    case 35  %Shekel7
        a=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
        c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
        f=0;
        for i=1:7
            f=f-1/((L-a(i,:))*(L-a(i,:))'+c(i));
        end
        
    case 36  %Shekel10
        a=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
        c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
        f=0;
        for i=1:10
            f=f-1/((L-a(i,:))*(L-a(i,:))'+c(i));
        end
        
    case 37  %Perm
        b=0.5;
        f1=0;
        for i=1:dim
            f2=0;
            for j=1:dim
                f2=f2+(j^i+b)*((L(j)/j)^i-1);
            end
            f1=f1+f2^2;
        end
        f=f1;
        
    case 38  %PowerSum
        dim=length(L);
        b=[8,18,44,114];
        f1=0;
        for i=1:dim
            f2=0;
            for j=1:dim
                f2=f2+L(j)^i;
            end
            f1=f1+(f2-b(i))^2;
        end
        f=f1;
        
    case 39  %Hartman
        a=[3 10 30;.1 10 35;3 10 30;.1 10 35];
        c=[1 1.2 3 3.2];
        p=[0.3689 0.117 0.2673;0.4699 0.4387 0.747;0.1091 0.8732 0.5547;0.03815 0.5743 0.8828];
        f=0;
        for i=1:4
            f=f-c(i)*exp(-(sum(a(i,:).*((L-p(i,:)).^2))));
        end
        
    case 40  %Hartman6
        a=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
        c=[1 1.2 3 3.2];
        p=[0.1312 0.1696 0.5569 0.0124 0.8283 0.5886;0.2329 0.4135 0.8307 0.3736 0.1004 0.9991;...
            0.2348 0.1415 0.3522 0.2883 0.3047 0.6650;0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];
        f=0;
        for i=1:4
            f=f-c(i)*exp(-(sum(a(i,:).*((L-p(i,:)).^2))));
        end
        
    case 41  %Griewank
        dim=length(L);
        f=sum(L.^2)/4000-prod(cos(L./sqrt([1:dim])))+1;
        
    case 42  %Ackley
        dim=length(L);
        f=-20*exp(-.2*sqrt(sum(L.^2)/dim))-exp(sum(cos(2*pi.*L))/dim)+20+exp(1);
        
    case 43  %Penalized
        dim=length(L);
        f=(pi/dim)*(10*((sin(pi*(1+(L(1)+1)/4)))^2)+sum((((L(1:dim-1)+1)./4).^2).*...
            (1+10.*((sin(pi.*(1+(L(2:dim)+1)./4)))).^2))+((L(dim)+1)/4)^2)+sum(100.*...
            ((L-10).^4).*(L>10)+100.*((-L-10).^4).*(L<(-10)));
        
    case 44  %Penalized2
        dim=length(L);
        f=0.1*((sin(3*pi*L(1)))^2+sum((L(1:dim-1)-1).^2.*(1+(sin(3.*pi.*L(2:dim))).^2))+...
            ((L(dim)-1)^2)*(1+(sin(2*pi*L(dim)))^2))+sum(100.*...
            ((L-10).^4).*(L>10)+100.*((-L-10).^4).*(L<(-10)));
        
    case 45  %Langerman2
        dim=length(L);
        para=parameter;
        a=para(:,1:10);
        c=para(:,11);
        f1=0;
        for i=1:dim
            f0=0;
            for j=1:dim
                f0=f0+(L(j)-a(i,j))^2;
            end
            s1=exp(-1/pi*f0);
            s2=cos(pi*f0);
            f1=f1-c(i)*s1*s2;
        end
        f=f1;
        
    case 46   %Langerman5
        dim=length(L);
        para=parameter;
        a=para(:,1:10);
        c=para(:,11);
        f1=0;
        for i=1:dim
            f0=0;
            for j=1:dim
                f0=f0+(L(j)-a(i,j))^2;
            end
            s1=exp(-1/pi*f0);
            s2=cos(pi*f0);
            f1=f1-c(i)*s1*s2;
        end
        f=f1;
        
    case 47  %Langerman10
        dim=length(L);
        para=parameter;
        a=para(:,1:10);
        c=para(:,11);
        f1=0;
        for i=1:dim
            f0= 0;
            for j=1:dim
                f0=f0+(L(j)-a(i,j))^2;
            end
            s1=exp(-1/pi*f0);
            s2=cos(pi*f0);
            f1=f1-c(i)*s1*s2;
        end
        f=f1;
        
    case 48  %FletcherPowell2
        [a,b,alpha]=a_b_alpha;
        dim=length(L);
        f=0;
        for i=1:dim
            A=0;
            B=0;
            for j=1:dim
                A=A+(a(i,j)*sin(alpha(j))+b(i,j)*cos(alpha(j)));
                B=B+(a(i,j)*sin(L(j))+b(i,j)*cos(L(j)));
            end
            f=f+(A-B)^2;
        end
        
    case 49   %FletcherPowell5
        [a,b,alpha]=a_b_alpha;
        dim=length(L);
        f=0;
        for i=1:dim
            A=0;
            B=0;
            for j=1:dim
                A=A+(a(i,j)*sin(alpha(j))+b(i,j)*cos(alpha(j)));
                B=B+(a(i,j)*sin(L(j))+b(i,j)*cos(L(j)));
            end
            f=f+(A-B)^2;
        end
        
    otherwise  %FletcherPowell10
        [a,b,alpha]=a_b_alpha;
        dim=length(L);
        f=0;
        for i=1:dim
            A=0;
            B=0;
            for j=1:dim
                A=A+(a(i,j)*sin(alpha(j))+b(i,j)*cos(alpha(j)));
                B=B+(a(i,j)*sin(L(j))+b(i,j)*cos(L(j)));
            end
            f=f+(A-B)^2;
        end
end
