function [edf] = idEDF(termo,edfs)
% idEDF Identifica se algum termo do vetor edfs existe na expressão
% simbolica termo.
% Entrada
% termo: Termo simbólico 
% edfs:  Vetor com o "nome" simbólico da edf a ser identificada no termo
% Saída
% edf: Índice do vetor edfs da edf identificada no termo, caso não encontre
% nenhuma retorna 0
edf = 0;
for k = 1:length(edfs)
    if has(termo,edfs(k))
        edf = k;
        break;
    end
end
end

