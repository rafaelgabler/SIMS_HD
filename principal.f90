!********************************************************************************************************!
!		            SUBROTINA PRINCIPAL DO PROGRAMA SIMS				         !
!											  		 !
! Ultima atualizacao em 24/06/2022							    		 !
! 											    		 !
!											    		 !
! Escrito por: Rafael Gabler Gontijo, PhD 						    		 !
!********************************************************************************************************!


! Esta subrotina tem como objetivo resolver as equacoes do momento linear e angular para uma
! suspensao de particulas rigidas monodispersas, na qual as seguintes forcas estao presentes:

! - Forcas de repulsao entre particulas (lubrificacao);
! - Forcas de contato entre particulas (Hertz);
! - Forcas por interacao magnetica entre os momentos de dipolo das particulas;
! - Forcas por interacao hidrodinamicas - sistemas periodicos

! Para a equacao do momento angular, consideram-se os seguintes torques presentes:

! - Torques magneticos (por interacao entre particula e por campo externo);
! - Torque viscoso;

! As escalas caracteristicas para o processo de adimensionalizacao sao:

! - Velocidade de Stokes
! - Tempo (a/U_s) ou (a²/D)

! As forças e torques por interacoes magneticas podem ser computadas de forma periodica ou nao


subroutine principal

use variaveis
use funcoes

! Numero de passos de tempo e passo de tempo

pi = acos(-1.0)
dt_inicial=min(0.01, pi/freqcampo)


! Se a simulacao e estatica (Monte Carlo) o numero de passos de tempo e 2, caso contrario depende do
! tempo total de simulacao escolhido pelo usuario

if(estatica) then
npast=2
else
npast=tempo/dt_inicial
end if

! Quantidade de numeros randomicos necessarios em cada timestep

nnr=3*N*rea

!********************************************************************************************************!

! Alocando variaveis na memoria

allocate(dx2(rea))
allocate(temposed(rea))
allocate(contador(rea))
allocate(posinicial(rea,2))
allocate(posfinal(rea,2))
allocate(X(rea,N,3))
allocate(XANTES(rea,N,3))
allocate(U(rea,N,3))
allocate(W(rea,N,3))
if(browniano)then
allocate(FORCAS(6,rea,N,3))
allocate(TORQUES(3,rea,N,3))
else
allocate(FORCAS(5,rea,N,3))
allocate(TORQUES(2,rea,N,3))
end if
allocate(FT(rea,N,3))
allocate(Tt(rea,N,3))
allocate(Di(rea,N,3))
allocate(aux1(rea,N))
allocate(aux2(rea,N))
allocate(aux3(rea,N))
allocate(aux4(rea,N))
allocate(nr(nnr))
allocate(dt(rea,N))
allocate(hidrodinamica_aux1(N,3))
allocate(hidrodinamica_aux2(N,3))
allocate(contribuicao_self(rea,N))
allocate(contribuicao_fisico(rea,N))
allocate(contribuicao_reciproco(rea,N))
allocate(hidro1(nb,3))
allocate(hidro2(nbr,3))
if(tmagper)then
allocate(auxt(N,3))
allocate(torquereal(nb,3))
allocate(torquereciproco(nbr,3))
allocate(cof4(2,10000))
allocate(cof5(2,10000))
allocate(cof7(2,10000))
end if
if(fmagper) then
allocate(cof6(2,10000))
allocate(cof8(2,10000))
allocate(auxf(N,3))
allocate(forcareal(nb,3))
allocate(forcareciproca(nbr,3))
end if
allocate(ILF(nb,3))
allocate(ILR(nbr,3))
allocate(XI(nb,rea,N,3))
allocate(cof1(2,10000))
allocate(cof2(2,10000))
allocate(cof3(2,10000))
if(leito)then
allocate(usistema(N,3))
end if
if(grafmag)then
allocate(magtempo(3,npast))
allocate(flutmag(N,rea))
end if
allocate(tempototal(npast))
if(agregado_inicial) then
allocate(centro_massa(rea,3))
end if
allocate(DIAM(rea,N))
allocate(beta(rea,N))
allocate(diarand(rea*N))
if(greenkubo)then
allocate(potencial(rea,N))
allocate(energia(rea,N))
allocate(energiaantes(rea,N))
allocate(heatcurrent(rea,N,3))
allocate(dedt(rea,N))
end if

!********************************************************************************************************!
!********************************************** ZERANDO TUDO ********************************************!
!********************************************************************************************************!

509 FORMAT(F30.4,F30.4,F30.4,F30.4)
666 FORMAT(F30.4,F30.4,F30.4,F30.4,F30.4,F30.4,F30.4)
1012 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)

X=0.0
U=0.0
aux1=0.0
aux2=0.0
aux3=0.0
Di=0.0
nr=0.0
FORCAS=0.0
FT=0.0
hidrodinamica_aux1=0.0
hidrodinamica_aux2=0.0
hidro1=0.0
hidro2=0.0
contribuicao_self=0.0
contribuicao_fisico=0.0
contribuicao_reciproco=0.0
torquereal=0.0
torquereciproco=0.0
auxt=0.0
shearratei=shearrate
contador=0

! Definindo o numero de Peclet rotacional em termos do translacional utilizando a relacao entre os coe-
! ficientes de difusao de Stokes-Einstein translacionais e rotacionais

Per=(4.0/3.0)*Pe

! Definindo o passo de tempo de todas as particulas igual ao inicial

do j=1,rea
do i=1,N
dt(j,i)=dt_inicial
end do
end do

! Definindo o tamanho das particulas

if(polidispersidade)then

call randomica(1.5,2.5,diarand,(N*rea),8)

do j=1,rea
do i=1,N
DIAM(j,i)=diarand((j-1)*N + i)
end do
end do

else

do j=1,rea
do i=1,N
DIAM(j,i)=2.0
end do
end do

end if

do j=1,rea
do i=1,n
beta(j,i) = DIAM(j,i)/DIAM(j,1)
end do
end do


! Calculando o tamanho do box para atender a fracao volumetrica desejada

if(agregado_inicial) then

ragreg=(N/phi)**(1.0/3.0)
l=10.0*ragreg
h=l

else
 l=3.0*(((N-2)/(razao*phi))*(4.0*3.1416)/(3.0))**(1.0/3.0)
 razao2=razao
 h=razao2*l
end if


! Criacao dos arquivos necessarios para armazenamento de resultados

if(.not.continua) then
 call gera_arquivos(posicao,velocidade,rea)
end if

! Calculando o numero pi

 pi=acos(-1.0)

! Definindo o parametro de convergencia da soma periodica qsi

qsi=1.0*((pi**0.5)/((l*l*h)**(1.0/3.0)))

! Tabelando as funcoes utilizadas na computacao das interacoes periodicas e criando os indices das 
! Lattices periodicas (caso algum tipo de interacao considere sistemas periodicos)

if(periodicidade) then
 call tabelagreen(qsi,l,nb,nbr,h)
 call estrutura_periodica
end if

! Distribuindo agora os momentos de dipolo iniciais das particulas

if(.not.continua) then
 call distribui_dipolo(Di,rea,N)
end if

! Vamos agora realizar uma simulacao paralela para todas as realizacoes

! Distribuindo inicialmente as particulas

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                         GENERATING INITIAL CONDITIONS                      *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

 call condicao_inicial


! Montando o diagrama de reversibilidade inicial

write(7*rea,*) 'zone t="Inicial"'
do j=1,rea
write(7*rea,*) X(j,N,1), X(j,N,2)
posinicial(j,1)= X(j,N,1)
posinicial(j,2)= X(j,N,2)
end do

write(7*rea,*) 'zone t="Final"'
if(posicao)then
do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
end do
end do
end if

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                   INITIAL CONDITIONS SUCCESSFULLY GENERATED                *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''


!********************************************************************************************************!
!********************************* INICIO DO PROCESSO NUMERICO DE SIMULACAO *****************************!
!********************************************************************************************************!


if(continua)then
iter=iter
auxiliar_continua=npast
else
iter=1
auxiliar_continua=npast-1
end if


aux_real=auxiliar_continua

! Iniciando o processo global numerico

if(printphi)then
call campo_phi(rea,k)   
end if

! Aqui é o looping que faz o tempo evoluir

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                                SIMULATING                                  *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

do k=iter, auxiliar_continua

k_real=k

if(shear) then

if(oscillatory) then
shearrate=shearratei*sin((1.0/freq)*dt_inicial*k)
end if

end if

! Calculando a energia das partículas antes de evoluir os potenciais e velocidades das particulas

if(greenkubo) then
do j=1,rea
do i=1,N
energiaantes(j,i) = potencial(j,i) !+ 0.5*(U(j,i,1)**2.0 + U(j,i,2)**2.0 + U(j,i,3)**2.0)**0.5
end do
end do
end if



! Primeira coisa de tudo: determinar o passo de tempo de cada particula antes de comecar a simulacao

do j=1,rea
do i=1,N
do q=1,N
if(i.ne.q) then
! Calcula-se a distancia entre a particula em questao e as outras particulas
r=(((X(j,i,1)-X(j,q,1))**2.0)+((X(j,i,2)-X(j,q,2))**2.0)+((X(j,i,3)-X(j,q,3))**2.0))**0.5
dt(j,i)=dt_inicial
dt(j,q)=dt_inicial
end if
end do
end do
end do

!********************************************************************************************************!
!********************************* Determinando as forcas Brownianas  ***********************************!
!********************************************************************************************************!

if(browniano)then
 call brownian
end if

!********************************************************************************************************!
!****************** Determinando as forcas de repulsao entre as particulas do box  **********************!
!********************************************************************************************************!

 call repulsao

!********************************************************************************************************!
!******************** Determinando as forcas de contato entre as particulas do box  *********************!
!********************************************************************************************************!

 call contato

!********************************************************************************************************!
!************** Determinando as forcas por interacao magnetica entre as particulas do box ***************!
!********************************************************************************************************!

if(.not.fmagper) then
 call forca_magnetica
end if

!********************************************************************************************************!
!************** Determinando as forcas por interacao magnetica entre as particulas e campo***************!
!********************************************************************************************************!

!if(externo) then
! call campo_externo
!else
do j=1,rea
do i=1,N
FORCAS(5,j,i,1)=0.0
FORCAS(5,j,i,2)=0.0
FORCAS(5,j,i,3)=0.0
end do
end do
!end if

! Impondo gravidade ou condicao de particulas neutrally buoyant

if(gravidade)then
do j=1,rea
do i=1,N
FORCAS(3,j,i,1)=0.0
FORCAS(3,j,i,2)=0.0
FORCAS(3,j,i,3)=-beta(j,i)**3.0
end do
end do
else
do j=1,rea
do i=1,N
FORCAS(3,j,i,1)=0.0
FORCAS(3,j,i,2)=0.0
FORCAS(3,j,i,3)=0.0
end do
end do
end if

! Implementação - Yureba
do j=1,rea
FORCAS(3,j,N,1)=0.0
FORCAS(3,j,N,2)=0.0
FORCAS(3,j,N,3)=-densityratio

FORCAS(3,j,N-1,1)=0.0
FORCAS(3,j,N-1,2)=0.0
FORCAS(3,j,N-1,3)=-1.0
end do


! Calculando todas as forcas que atuam nas particulas

do j=1,rea
do i=1,N
FT(j,i,1)=FORCAS(1,j,i,1)+FORCAS(2,j,i,1)+FORCAS(3,j,i,1)+FORCAS(4,j,i,1)+FORCAS(5,j,i,1)+FORCAS(6,j,i,1)
FT(j,i,2)=FORCAS(1,j,i,2)+FORCAS(2,j,i,2)+FORCAS(3,j,i,2)+FORCAS(4,j,i,2)+FORCAS(5,j,i,2)+FORCAS(6,j,i,2)
FT(j,i,3)=FORCAS(1,j,i,3)+FORCAS(2,j,i,3)+FORCAS(3,j,i,3)+FORCAS(4,j,i,3)+FORCAS(5,j,i,3)+FORCAS(6,j,i,3)
end do
end do

!********************************************************************************************************!
!******************************* DETERMINANDO AS INTERACOES PERIODICAS  *********************************!
!********************************************************************************************************!

if(periodicidade) then
call intper
end if

if(.not.ligaih) then
if(inertia) then

do j=1,rea
do i=1,N
 call resvel(U(j,i,1),dt_inicial,St,FT(j,i,1))
 call resvel(U(j,i,2),dt_inicial,St,FT(j,i,2))
 call resvel(U(j,i,3),dt_inicial,St,FT(j,i,3))
end do
end do
else
U=FT
endif
end if

! Determinando a velocidade em caso de cisalhamento simples

if(shear)then
do j=1,rea
do i=1,N
U(j,i,2) = U(j,i,2) + shearrate*X(j,i,3)
end do
end do
end if


!***************************** MODELO DE LEITO FLUIDIZADO MAGNETICO ***********************************!

if(leito)then
do i=1,N
usistema(i,1)=sum(U(:,i,1))/rea
usistema(i,2)=sum(U(:,i,2))/rea
usistema(i,3)=sum(U(:,i,3))/rea
end do


do q=1,rea
do i=1,N
U(q,i,1)=U(q,i,1)-usistema(i,1)
U(q,i,2)=U(q,i,2)-usistema(i,2)
U(q,i,3)=U(q,i,3)-usistema(i,3)
end do
end do
end if

! Calculando a energia total (cinética + potencial) das partículas para uso em Green-Kubo 

if(greenkubo) then
do j=1,rea
do i=1,N
energia(j,i) = potencial(j,i) !+ 0.5*(U(j,i,1)**2.0 + U(j,i,2)**2.0 + U(j,i,3)**2.0)**0.5

dedt(j,i)= (energiaantes(j,i)-energia(j,i))/dt_inicial

heatcurrent(j,i,1) = energia(j,i)*U(j,i,1) + X(j,i,1)*dedt(j,i)
heatcurrent(j,i,2) = energia(j,i)*U(j,i,2) + X(j,i,2)*dedt(j,i)
heatcurrent(j,i,3) = energia(j,i)*U(j,i,3) + X(j,i,3)*dedt(j,i)
end do
end do
end if


!******************************************************************************************************!

! Calculando a posicao atual das particulas (Utilizando Euler)

XANTES=X

do j=1,rea
 do i=1,N
 call respos(X(j,i,1),dt(j,i),U(j,i,1))
 call respos(X(j,i,2),dt(j,i),U(j,i,2))
 call respos(X(j,i,3),dt(j,i),U(j,i,3))
 end do
end do


deslocamentox = sum(X(:,:,1) - XANTES(:,:,1))/(N*rea)
deslocamentoy = sum(X(:,:,2) - XANTES(:,:,2))/(N*rea)
deslocamentoz = sum(X(:,:,3) - XANTES(:,:,3))/(N*rea)

deslocamentomedio = sum(((X(:,:,1) - XANTES(:,:,1))**2.0 + (X(:,:,2) - XANTES(:,:,2))**2.0 + (X(:,:,3) - XANTES(:,:,3))**2.0)**0.5)/(N*rea)

write(300*rea,*) deslocamentox, deslocamentoy, deslocamentoz, deslocamentomedio, k*dt_inicial

!******************************************************************************************************!
! Implementando a condicao de contorno de periodicidade

!if(periodicidade) then
do j=1,rea
do i=1,N
if(X(j,i,1).gt.l) then
X(j,i,1)=X(j,i,1)-l
end if
if(X(j,i,1).lt.0.0)then
X(j,i,1)=l-X(j,i,1)
end if
if(X(j,i,2).gt.l) then
X(j,i,2)=X(j,i,2)-l
end if
if(X(j,i,2).lt.0.0)then
X(j,i,2)=l-X(j,i,2)
end if
if(X(j,i,3).gt.h) then
X(j,i,3)=X(j,i,3)-h
end if
if(X(j,i,3).lt.0.0)then
X(j,i,3)=h-X(j,i,3)
end if
end do
end do
!end if

!********************************************************************************************************!

! Montando o diagrama de reversibilidade final

do j=1,rea
if(X(j,N,3).le.((h-laux-laux)/2.0))then
if(j.ne.contador(j)) then
write(7*rea,*) X(j,N,1), X(j,N,2)
contador(j)=j
temposed(j)=k*dt_inicial
posfinal(j,1)= X(j,N,1)
posfinal(j,2)= X(j,N,2)
end if
end if
end do

! Escrevendo em um arquivo de saida a posicao e velocidade de cada particula em cada instante de tempo


if(continua)then

teste1=k/n3
teste2=k/n2

if(teste1.eq.teste2) then


if(posicao)then
do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
end do
end do
end if



write(100*rea,'(A12,I6,A1)') 'zone t="',k,'"'
do j=1,rea
write(100*rea,*) X(j,N,1),X(j,N,2),X(j,N,3), Di(j,N,1), Di(j,N,2), Di(j,N,3)
end do

write(*,*) 'Percentual:', (k_real/aux_real)*100, '%'

if(agregado_inicial) then
write(666,*) k*dt_inicial, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if


if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if


if(greenkubo)then
do j=1,rea
write(2012*j,*) 'zone t="',k,'"'
do i=1,N
write(2012*j,1012) potencial(j,i),energia(j,i),dedt(j,i),heatcurrent(j,i,1),heatcurrent(j,i,2),heatcurrent(j,i,3) 
end do
end do
end if


end if



else

if(k.eq.1) then

if(posicao)then
do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
if(gravadipolo)then
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
else
write(j,509)X(j,i,1),X(j,i,2),X(j,i,3),DIAM(j,i)
end if
end do
end do
end if

write(100*rea,'(A12,I6,A1)') 'zone t="',k,'"'
do j=1,rea
write(100*rea,*) X(j,N,1),X(j,N,2),X(j,N,3), Di(j,N,1), Di(j,N,2), Di(j,N,3)
end do

write(*,*) 'Percentual:', (k_real/aux_real)*100, '%'

if(agregado_inicial) then
write(666,*) k*dt_inicial, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if


if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if

end if

teste1=k/n3
teste2=k/n2

if(k.ne.1) then

if(teste1.eq.teste2) then
if(posicao)then
do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
if(gravadipolo)then
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
else
write(j,509)X(j,i,1),X(j,i,2),X(j,i,3),DIAM(j,i)
end if
end do
end do
end if

write(100*rea,'(A12,I6,A1)') 'zone t="',k,'"'
do j=1,rea
write(100*rea,*) X(j,N,1),X(j,N,2),X(j,N,3), Di(j,N,1), Di(j,N,2), Di(j,N,3)
end do

write(*,*) 'Percentual:', (k_real/aux_real)*100, '%'  

if(agregado_inicial) then
write(666,*) k*dt_inicial, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if

if(greenkubo)then
do j=1,rea
write(2012*j,*) 'zone t="',k,'"'
do i=1,N
write(2012*j,1012) potencial(j,i),energia(j,i),dedt(j,i),heatcurrent(j,i,1),heatcurrent(j,i,2),heatcurrent(j,i,3)
end do
end do
end if



if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if
end if

end if
end if

 
!********************************************************************************************************!
!********************************************************************************************************!
!***************************** SOLUCAO DO MOVIMENTO ROTACIONAL DAS PARTICULAS ***************************!
!********************************************************************************************************!
!********************************************************************************************************!

if(torque) then

!*******************************************************************************************!
!**************** Determinando os torques magneticos sobre as particulas *******************!
!*******************************************************************************************!

 if(.not.tmagper)then
 call torque_magnetico
 end if

if(externo) then

if(rotating)then

call rotating_field(alpha, freqcampo*k*dt_inicial)

else

if(oscilacampo) then
 call torque_externo(alpha*sin(freqcampo*k*dt_inicial))
else
 call torque_externo(alpha)
end if
end if

end if

! Caso existam torques brownianos eles ja foram computados acima, na subrotina "brownian"

! Computando os torques totais que atuam sobre as particulas

do j=1,rea
do i=1,N
if(browniano) then
Tt(j,i,1)= TORQUES(1,j,i,1) + TORQUES(2,j,i,1) + TORQUES(3,j,i,1)
Tt(j,i,2)= TORQUES(1,j,i,2) + TORQUES(2,j,i,2) + TORQUES(3,j,i,2)
Tt(j,i,3)= TORQUES(1,j,i,3) + TORQUES(2,j,i,3) + TORQUES(3,j,i,3)
else
Tt(j,i,1)= TORQUES(1,j,i,1) + TORQUES(2,j,i,1) 
Tt(j,i,2)= TORQUES(1,j,i,2) + TORQUES(2,j,i,2) 
Tt(j,i,3)= TORQUES(1,j,i,3) + TORQUES(2,j,i,3) 
end if
end do
end do


! Resolvendo a velocidade angular

if(mistura)then
do j=1,rea
do i=(percentual*N)+1,N
 call resomega(W(j,i,1),dt(j,i),Str,Tt(j,i,1))
 call resomega(W(j,i,2),dt(j,i),Str,Tt(j,i,2))
 call resomega(W(j,i,3),dt(j,i),Str,Tt(j,i,3))
end do
end do
else
do j=1,rea
do i=1,N
 call resomega(W(j,i,1),dt(j,i),Str,Tt(j,i,1))
 call resomega(W(j,i,2),dt(j,i),Str,Tt(j,i,2))
 call resomega(W(j,i,3),dt(j,i),Str,Tt(j,i,3))
! call resomega_sem_inercia(W(j,i,1),Tt(j,i,1))
! call resomega_sem_inercia(W(j,i,2),Tt(j,i,2))
! call resomega_sem_inercia(W(j,i,3),Tt(j,i,3))
end do
end do
end if

if(shear)then
do j=1,rea
do i=1,N
W(j,i,1) = W(j,i,1) - shearrate*0.5
end do
end do
end if

! Evoluindo o vetor momento de dipolo magneticos das particulas

if(mistura)then
do j=1,rea
do i=(N*percentual)+1,N
 call evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,3),W(j,i,2),W(j,i,3),dt(j,i))
 call evoldip(Di(j,i,2),Di(j,i,3),Di(j,i,1),W(j,i,3),W(j,i,1),dt(j,i))
 call evoldip(Di(j,i,3),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt(j,i))
end do
end do
else
do j=1,rea
do i=1,N
 call evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,3),W(j,i,2),W(j,i,3),dt(j,i))
 call evoldip(Di(j,i,2),Di(j,i,3),Di(j,i,1),W(j,i,3),W(j,i,1),dt(j,i))
 call evoldip(Di(j,i,3),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt(j,i))
end do
end do
end if

! Normalizando os dipolos

do j=1,rea
if(mistura)then
do i=1,(percentual*N)
Di(j,i,1)=0.0
Di(j,i,2)=0.0
Di(j,i,3)=0.0
end do
do i=(percentual*N)+1,N
Di(j,i,1)=Di(j,i,1)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
Di(j,i,2)=Di(j,i,2)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
Di(j,i,3)=Di(j,i,3)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
end do
else
do i=1,N
Di(j,i,1)=Di(j,i,1)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
Di(j,i,2)=Di(j,i,2)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
Di(j,i,3)=Di(j,i,3)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
end do
end if
end do

! Determinando a magnetizacao da suspensao, caso seja solicitado pelo usuario

if(grafmag) then
 call media_ativa(Di,N,rea,magtempo(1,k),1)
 call media_ativa(Di,N,rea,magtempo(2,k),2)
 call media_ativa(Di,N,rea,magtempo(3,k),3)
! Calculando a barra de erro da magnetizacao de equilibrio
do j=1,rea
do i=1,N
flutmag(i,j)=(Di(j,i,3)-magtempo(3,k))**2.0
end do
end do
erromag=((1.0/(N*rea))*sum(flutmag))**0.5
!write(5*rea,*)magtempo(3,k),k*dt_inicial, ((magtempo(k)-magtempo(k-1))/dt_inicial),  erromag
write(5*rea,*)magtempo(1,k),magtempo(2,k),magtempo(3,k), (magtempo(1,k)**2.0+magtempo(2,k)**2.0+magtempo(3,k)**2.0 )**0.5, k*dt_inicial
end if

tempototal(k)=k*dt_inicial
end if
end do

! Analisando aqui os deslocamentos quadraticos medios por meio das diferenças nos diagramas 
! de reversibilidade (inicial e final)

do j=1,rea
dx2(j)= (posinicial(j,1)-posfinal(j,1))**2.0 + (posinicial(j,2)-posfinal(j,2))**2.0
end do

! Calculando o deslocamento medio quadratico da particula N

desquadrado=sum(dx2)/rea
tempomedio=sum(temposed)

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                         PARTICLE DISPERSION REPORT	                      *'
print *,'*                                                                            *'
print *,'******************************************************************************'
write(*,*)''
write(*,*)'MEAN SQUARE DISPLACEMENT OF THE UPPER PARTICLE:',desquadrado
write(*,*)'MEAN SEDIMENTATION TIME OF THE UPPER PARTICLE:',tempomedio
write(*,*)'PARTICLE HYDRODYNAMIC DIFFUSION COEFFICIENT:', (desquadrado/tempomedio)

write(23802,*)'MEAN SQUARE DISPLACEMENT OF THE UPPER PARTICLE:',desquadrado
write(23802,*)'MEAN SEDIMENTATION TIME OF THE UPPER PARTICLE:',tempomedio
write(23802,*)'PARTICLE HYDRODYNAMIC DIFFUSION COEFFICIENT:', (desquadrado/tempomedio)



if(printphi) then
call campo_phi(rea,k)
end if



if(fator)then
 call fator_estrutura(X,N,l,h,dt_inicial,rea)
end if



! Fechando os arquivos abertos
do j=1,2*rea
 close(j)
end do

 close(100*rea)
 close(300*rea)

if(greenkubo) then
do j=1,rea
 close(2012*j)
end do
end if

! Desalocando todas as matrizes e vetores

deallocate(X, STAT = DeAllocateStatus)
deallocate(U, STAT = DeAllocateStatus)
deallocate(FORCAS, STAT = DeAllocateStatus)
deallocate(FT, STAT = DeAllocateStatus)
deallocate(nr, STAT = DeAllocateStatus)
deallocate(dt, STAT = DeAllocateStatus)
deallocate(hidrodinamica_aux1, STAT = DeAllocateStatus)
deallocate(hidrodinamica_aux2, STAT = DeAllocateStatus)
deallocate(hidro1, STAT = DeAllocateStatus)
deallocate(hidro2, STAT = DeAllocateStatus)
deallocate(ILF, STAT = DeAllocateStatus)
deallocate(ILR, STAT = DeAllocateStatus)
deallocate(XI, STAT = DeAllocateStatus)
deallocate(Tt, STAT = DeAllocateStatus)
deallocate(Di, STAT = DeAllocateStatus)
deallocate(aux1, STAT = DeAllocateStatus)
deallocate(aux2, STAT = DeAllocateStatus)
deallocate(aux3, STAT = DeAllocateStatus)
deallocate(aux4, STAT = DeAllocateStatus)
deallocate(contribuicao_self, STAT = DeAllocateStatus)
deallocate(contribuicao_fisico, STAT = DeAllocateStatus)
deallocate(contribuicao_reciproco, STAT = DeAllocateStatus)
if(tmagper)then
deallocate(auxt, STAT = DeAllocateStatus)
deallocate(torquereal, STAT = DeAllocateStatus)
deallocate(torquereciproco, STAT = DeAllocateStatus)
deallocate(cof4, STAT = DeAllocateStatus)
deallocate(cof5, STAT = DeAllocateStatus)
deallocate(cof7, STAT = DeAllocateStatus)
end if
if(fmagper) then
deallocate(cof6, STAT = DeAllocateStatus)
deallocate(cof8, STAT = DeAllocateStatus)
deallocate(auxf, STAT = DeAllocateStatus)
deallocate(forcareal, STAT = DeAllocateStatus)
deallocate(forcareciproca, STAT = DeAllocateStatus)
end if
deallocate(ILF, STAT = DeAllocateStatus)
deallocate(ILR, STAT = DeAllocateStatus)
deallocate(XI, STAT = DeAllocateStatus)
deallocate(cof1, STAT = DeAllocateStatus)
deallocate(cof2, STAT = DeAllocateStatus)
deallocate(cof3, STAT = DeAllocateStatus)
if(leito)then
deallocate(usistema, STAT = DeAllocateStatus)
end if
if(grafmag)then
deallocate(magtempo, STAT = DeAllocateStatus)
deallocate(flutmag, STAT = DeAllocateStatus)
end if
deallocate(tempototal, STAT = DeAllocateStatus)
if(agregado_inicial) then
deallocate(centro_massa, STAT = DeAllocateStatus)
end if
deallocate(DIAM, STAT = DeAllocateStatus)
deallocate(beta, STAT = DeAllocateStatus)
deallocate(diarand, STAT = DeAllocateStatus)
if(greenkubo)then
deallocate(potencial, STAT = DeAllocateStatus)
deallocate(energia, STAT = DeAllocateStatus)
deallocate(energiaantes, STAT = DeAllocateStatus)
deallocate(heatcurrent, STAT = DeAllocateStatus)
deallocate(dedt, STAT = DeAllocateStatus)
end if







write(*,*) ''

write(*,*) 'Finalizando o modulo de processamento...'

write(*,*) ''

end subroutine principal
