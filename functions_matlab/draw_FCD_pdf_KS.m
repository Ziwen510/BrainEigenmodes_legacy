function ksStat = draw_FCD_pdf_KS(mat, version, name)
% compareUpperTriPdf Extracts upper triangular values from input matrix,
% computes its PDF over [0,1] with 100 bins, loads another PDF from .mat,
% plots both PDFs, computes KS statistic, annotates plot, and saves figure as PNG.

% Define bin edges and centers
binEdges = linspace(0,1,101);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

% Extract upper triangular values (excluding diagonal)
utValues = mat(triu(true(size(mat)),1));

% Compute PDF for input matrix
pdf1 = histcounts(utValues, binEdges, 'Normalization', 'pdf');

% Load the other PDF
otherPdfMatFile = sprintf("./data/empirical/group_FCD_pdf_%s.mat", version);
S = load(otherPdfMatFile);
if isfield(S,'pdf')
    pdf2 = S.pdf;
else
    fn = fieldnames(S);
    pdf2 = S.(fn{2});
end
% disp(pdf2);

if numel(pdf2) ~= numel(pdf1)
    error('Loaded PDF length (%d) does not match expected number of bins (%d).', numel(pdf2), numel(pdf1));
end


% Compute CDFs
binWidth = binEdges(2) - binEdges(1);
pdf2 = pdf2 / binWidth;
cdf1 = cumsum(pdf1) * binWidth;
cdf2 = cumsum(pdf2) * binWidth;

% Compute KS statistic\ nksStat = max(abs(cdf1 - cdf2));
ksStat = max(abs(cdf1 - cdf2));

% Plot PDFs
figure;
plot(binCenters, pdf1, 'LineWidth', 2); hold on;
plot(binCenters, pdf2, 'LineWidth', 2);
xlabel('Value');
ylabel('Density');
legend('Simulated PDF', 'Empirical PDF');
title(sprintf('%s with %s', name, version), 'Interpreter','none');

% Annotate KS statistic
xPos = min(binCenters) + 0.05*(max(binCenters) - min(binCenters));
yPos = max([pdf1, pdf2]) * 0.9;
text(xPos, yPos, sprintf('KS = %.3f', ksStat), 'FontSize', 12);

% Save figure as PNG
saveas(gcf, sprintf("./outputs/%s_FCD_pdf(%s).png", name, version));
end
