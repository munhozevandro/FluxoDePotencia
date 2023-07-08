#--------------------
# Fluxo potencia para sistema de distribuicao de energia eletrico
#--------------------

#--------------------
# Conjuntos
#--------------------
set Ob;					   # Conjunto de barras do sistema
set Ol within Ob cross Ob; # Conjunto de linhas do sistema

#--------------------
# Parametros de entrada do fluxo de potencia
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
var I{Ol}; # as correntes nas linhas
var P{Ol}; # potencia ativa nas linhas
var Q{Ol}; # potencia reativa nas linhas

var PS{Ob}; # potencia ativa da subestação
var QS{Ob}; # potencia reativa da subestação
var V{Ob};  # tensão nas barras


#--------------------
# Funcao objetivo
minimize FuncaoObjetivo: sum{ (i,j) in Ol}( R[i,j]*I[i,j]^2 );

#--------------------
#Conjunto de restrições
subject to BalancoPotenciaAtiva{i in Ob}:
	sum{(k,i) in Ol}(P[k,i]) - sum{(i,j) in Ol}(P[i,j] + R[i,j]*I[i,j]^2) + PS[i] = PD[i]; 

subject to BalancoPotenciaReativa{i in Ob}:
	sum{(k,i) in Ol}(Q[k,i]) - sum{(i,j) in Ol}(Q[i,j] + X[i,j]*I[i,j]^2) + QS[i] = QD[i];
	
subject to quedaTensao{(i,j) in Ol}:
	V[i]^2 -  2*(R[i,j]*P[i,j] + X[i,j]*Q[i,j]) - Z2[i,j]*I[i,j]^2 - V[j]^2 = 0;	
	
subject to potenciaAparente {(i,j) in Ol}:
	Vnom^2 * I[i,j]^2 = P[i,j]^2 + Q[i,j]^2;
	
subject to limiteCorrente{ (i,j) in Ol}:
	0 <= I[i,j] <=Imax[i,j];
	
subject to limiteTensao{i in Ob}:
	Vmin <= V[i] <= Vmax;

