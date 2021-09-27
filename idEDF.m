function [edf] = idEDF(termo,edfs)
% idEDF Identifica se algum termo do vetor edfs existe na express�o
% simbolica termo.
% Entrada
% termo: Termo simb�lico 
% edfs:  Vetor com o "nome" simb�lico da edf a ser identificada no termo
% Sa�da
% edf: �ndice do vetor edfs da edf identificada no termo, caso n�o encontre
% nenhuma retorna 0
edf = 0;
for k = 1:length(edfs)
    if has(termo,edfs(k))
        edf = k;
        break;
    end
end
end

