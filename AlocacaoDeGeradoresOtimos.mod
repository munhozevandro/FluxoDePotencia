#--------------------
# Fluxo potencia para sistema de distribuicao de energia eletrico
#--------------------

# Fluxo de Potencia - I²,V²
# Fluxo de Potencia Modificado - (I^sqrt, V^sqrt, Vnom² = Vj², mudando = para >=)
# Consequentemente faz-se:
# Planejamento de Bancos Capacitivos
# Planejamento deGeradores distribuidos



#--------------------
# Conjuntos
#--------------------

set Ob;					   # Conjunto de barras do sistema
set Ol within Ob cross Ob; # Conjunto de linhas do sistema

#--------------------
# Paramentros de entrada do fluxo de potencia
#--------------------
param Vnom;		#Tensão nominal da rede em KV
param Barra_SE;	# Barra da substeação
param Vmin;		# Tensão minima da rede
param Vmax; 	# Tensão máxima da rede

# Parametros das barras

param PD{Ob};	# potencia ativa nas barras
param QD{Ob};	# potencia reativa nas barras
param QC{Ob};	# capacitor caso exista na barra

# Parametros para as linhas

param linea{Ol};	#id da linha
param R{Ol};		# resistencia da linha
param X{Ol};		# reatancia da linha
param Imax{Ol}; 	# Corrente maxima da linha
param Z2{Ol};		# Z^2 = R^2 + X^2

#--------------------
# Variaveis
#--------------------
var Iquadrado{Ol}; # as correntes nas linhas ao quadrado
var P{Ol}; # potencia ativa nas linhas
var Q{Ol}; # potencia reativa nas linhas

var PS{Ob}; # potencia ativa da subestação
var QS{Ob}; # potencia reativa da subestação
var Vquadrado{Ob};  # tensão nas barras ao quadrado


# Para a alocacao de Geradores Otimos

var Pgd{Ob};
var Qgd{Ob};
var Wgd{Ob}, binary;

param Cinst_gd = 10000;
param Ngd_syst = 5;
param Pgd_max = 1000; #1000kW
param fpgd = 0.92;
param ke = 1000;

#--------------------
# Funcao objetivo
minimize FuncaoObjetivo: ke * sum{ (i,j) in Ol}( R[i,j]*Iquadrado[i,j] ) + sum{i in Ob}(Cinst_gd * Wgd[i]);

#--------------------
# Conjunto de restrições
subject to BalancoPotenciaAtiva{i in Ob}:
	sum{(k,i) in Ol}(P[k,i]) - sum{(i,j) in Ol}(P[i,j] + R[i,j]*Iquadrado[i,j]) + PS[i] + Pgd[i] = PD[i]; 

subject to BalancoPotenciaReativa{i in Ob}:
	sum{(k,i) in Ol}(Q[k,i]) - sum{(i,j) in Ol}(Q[i,j] + X[i,j]*Iquadrado[i,j]) + QS[i] + Qgd[i]= QD[i];
	
subject to quedaTensao{(i,j) in Ol}:
	Vquadrado[i] -  2*(R[i,j]*P[i,j] + X[i,j]*Q[i,j]) - Z2[i,j]*Iquadrado[i,j] - Vquadrado[j] = 0;	
	
subject to potenciaAparente {(i,j) in Ol}:
	Vquadrado[j] * Iquadrado[i,j] >= P[i,j]^2 + Q[i,j]^2;
	
subject to limiteCorrentes{ (i,j) in Ol}:
	0 <= Iquadrado[i,j] <=Imax[i,j]^2;
	
subject to limiteTensao{i in Ob}:
	Vmin^2 <= Vquadrado[i] <= Vmax^2;

#--------------------
# equações do planejamento de capacitores 

subject to gerador_0{i in Ob}:
	Qgd[i] = Pgd[i]*tan(acos(fpgd));

subject to gerador_1{i in Ob}:
	Pgd[i] <= Pgd_max*Wgd[i];
	
subject to gerador_3:
	sum{i in Ob}(Wgd[i]) <= Ngd_syst;		

