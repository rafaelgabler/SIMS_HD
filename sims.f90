program magpart

use variaveis

! Titulo do programa e apresentacao do mesmo

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                  SIMS - SIMULATION OF MAGNETIC SUSPENSIONS	              *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                      HYDRODYNAMIC DISPERSION MODULE                        *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                     PROF. RAFAEL GABLER GONTIJO, PhD                       *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                        IN DEVELOPMENT SINCE 2009                           *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                         LAST UPDATE: 24/06/2022                            *'
print *,'*                                                                            *'
print *,'******************************************************************************'

print *,''
print *,'Numerical simulation of magnetic suspensions of hard spheres'
print *,''
print *,'Techniques implemented in this code:'
print *,''
print *,'1 - Langevin Dynamics'
print *,'2 - Stokesian Dynamics'
print *,'3 - Ewald Summation Technique'
print *,''
!**************************************************************************************************************************************!
 call entrada

! Começa a contar o tempo de simulacao (ativa o cronometro)

 call cpu_time(ti)

! Chama a subrotina principal que executa de fato a simulacao

 call principal

if(estatistica) then
 call saida
end if

! Para o cronometro

 call cpu_time(tf)

! Calcula o tempo de processamento

 tpros=tf-ti

write(*,*) ''
write(*,*) 'O tempo total de processamento foi de:',tpros,'segundos'
write(*,*) ''


end program magpart
