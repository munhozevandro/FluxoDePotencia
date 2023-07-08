#--------------------
# Fluxo potencia para sistema de distribuicao de energia eletrico
#--------------------

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

#--------------------
# Para o planejamento de capacitores

var Qfx{Ob}; 	# potencia reativa do capacitores
var Nfx{Ob} >=0, integer;
var Wfx{Ob}, binary;

param Qfx_unit = 300; # potencia reativa base do capacitor [kVAr]
param Nfx_max_system = 5; # numero maximo de capacitores do sistema
param Nfx_max = 4; # numero maximo de capacitores instalados na barra i 

#--------------------
# Custo de Instalação

param Cfx_inst = 10000;
param Cfx_unit = 9000; 
param ke = 1000;

#--------------------
# Funcao objetivo
minimize FuncaoObjetivo: ke * sum{ (i,j) in Ol}( R[i,j]*Iquadrado[i,j]) + sum{i in Ob}(Cfx_inst * Wfx[i] + Cfx_unit * Nfx[i]);
#--------------------
# Conjunto de restrições
subject to BalancoPotenciaAtiva{i in Ob}:
	sum{(k,i) in Ol}(P[k,i]) - sum{(i,j) in Ol}(P[i,j] + R[i,j]*Iquadrado[i,j]) + PS[i] = PD[i]; 

subject to BalancoPotenciaReativa{i in Ob}:
	sum{(k,i) in Ol}(Q[k,i]) - sum{(i,j) in Ol}(Q[i,j] + X[i,j]*Iquadrado[i,j]) + QS[i] + Qfx[i]= QD[i];
	
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

subject to capacitor_1{i in Ob}:
	Qfx[i] = Nfx[i]*Qfx_unit;
	
subject to capacitor_2{i in Ob}:
	Nfx[i] <= Wfx[i] * Nfx_max;
	
subject to capacitor_3:
	sum{i in Ob}(Wfx[i]) <= Nfx_max_system;		

