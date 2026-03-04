// =============================================================================
// Nuclei segmentation, measurement and thumbnail cropping
// Figures: Upstream of image_feature_analysis.Rmd (produces per-nucleus TIFFs)
//
// Processes max-projected TIFFs: segments nuclei from DAPI (channel 1)
// using Otsu thresholding + watershed, measures all channels per ROI,
// and crops/saves multi-channel thumbnails for downstream R analysis.
// Run in Fiji/ImageJ.
// =============================================================================

// Define directories
inputDir = getDirectory("Select Input Directory");
outputDir = inputDir + "processed" + File.separator;
if (!File.exists(outputDir)) File.makeDirectory(outputDir);

// Create thumbnails folder
thumbs_folder = outputDir + "MSL2gp" + File.separator;
if (!File.exists(thumbs_folder)) File.makeDirectory(thumbs_folder);

// Get user input for dataset name
gene_of_interest = getString("Specify name of data to be analyzed:", "MSL2gp");

// Get list of TIFF files
fileList = getFileList(inputDir);

// Set measurements
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=3");

// Process each image
setBatchMode(false);
for (i = 0; i < fileList.length; i++) {
    filePath = inputDir + fileList[i];
    if (endsWith(filePath, ".tif") || endsWith(filePath, ".tiff")) {
        // Open original image
        open(filePath);
        originalTitle = getTitle();
        baseName = replace(originalTitle, ".tif", "");
        baseName = replace(baseName, ".tiff", "");

        print("Processing image: " + originalTitle);

        // Duplicate channel 1 (DAPI) for processing
        run("Duplicate...", "duplicate channels=1");
        rename("DAPI_Mask");

        // Process DAPI channel to create mask
        selectWindow("DAPI_Mask");
        run("Median...", "radius=3");
        run("Subtract Background...", "rolling=125");
        setAutoThreshold("Otsu dark");
        run("Convert to Mask");
        run("Watershed");

        // Save mask
        selectWindow("DAPI_Mask");
        saveAs("Tiff", outputDir + baseName + "_mask.tif");

        // Analyse particles: size 12-50 µm², circularity 0.6-1.0
        run("Analyze Particles...", "size=12-50 circularity=0.6-1.0 exclude add");
        roiManager("Save", outputDir + baseName + "_ROIs.zip");

        // Switch to the original image
        selectWindow(originalTitle);

        // Get number of channels
        Stack.getDimensions(width, height, channels, slices, frames);

        // Measure for each ROI in all channels
        roiCount = roiManager("count");
        for (roiIndex = 0; roiIndex < roiCount; roiIndex++) {
            roiManager("Select", roiIndex);
            roiName = Roi.getName();

            for (channel = 1; channel <= channels; channel++) {
                Stack.setChannel(channel);
                run("Measure");
                setResult("ROI", nResults-1, roiName);
                setResult("Channel", nResults-1, channel);
                setResult("Image", nResults-1, originalTitle);
            }

            // Crop and save thumbnail with all channels
            Roi.getBounds(x, y, width, height);
            WidthA = 110;
            HeightA = 110;
            xa = x - WidthA/2 + width/2;
            ya = y - HeightA/2 + height/2;
            xa = maxOf(0, minOf(getWidth()-WidthA, xa));
            ya = maxOf(0, minOf(getHeight()-HeightA, ya));

            makeRectangle(xa, ya, WidthA, HeightA);
            run("Duplicate...", "duplicate");
            saveAs("Tiff", thumbs_folder + gene_of_interest + "_" + roiName + "_" + baseName + ".tif");
            close();
        }

        // Save the results to a CSV file
        saveAs("Results", outputDir + baseName + "_measurements.csv");
        print("Saved measurements: " + outputDir + baseName + "_measurements.csv");

        // Reset ROI Manager and clear results
        roiManager("reset");
        run("Clear Results");

        // Close all open images
        close("*");
    }
}

print("Batch processing complete.");
