
pdfFile = "MMR_Schematic9_4points_20240725_1.pdf"; % String
images = PDFtoImg(pdfFile);

function images = PDFtoImg(pdfFile)
import org.apache.pdfbox.*
import java.io.*
filename = fullfile(pwd,pdfFile);
jFile = File(filename);
document = pdmodel.PDDocument.load(jFile);
pdfRenderer = rendering.PDFRenderer(document);
count = document.getNumberOfPages();
images = [];
dpi_num = 600;
for ii = 1:count
    bim = pdfRenderer.renderImageWithDPI(ii-1, dpi_num, rendering.ImageType.RGB);
    images = [images (filename + "-" +"Page" + ii + ".png")];
    tools.imageio.ImageIOUtil.writeImage(bim, filename + "-" +"Page" + ii + ".png", 300);
end
document.close()

end