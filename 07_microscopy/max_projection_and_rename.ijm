// =============================================================================
// Max-intensity Z-projection and batch renaming
// Figures: Upstream of image_feature_analysis.Rmd (produces single-plane TIFFs)
//
// Opens .lif stacks, performs max-intensity Z-projection, renames with
// batch number and saves as TIFF.
// Run in Fiji/ImageJ with open .lif images.
// =============================================================================

// Prompt user for the name of the subfolder to save processed images
subfolder = getString("Enter subfolder name for processed images:", "max_projection_");

// Select the parent folder; the subfolder will be created within the parent folder
outputFolder = getDirectory("Choose directory") + subfolder + File.separator;
File.makeDirectory(outputFolder);
print("Output folder: " + outputFolder);

// Prompt for batch number
batchNumber = getNumber("Enter the batch number:", 1);

// Get list of all open images (these are the images to process)
imageList = getList("image.titles");

// Function to generate a timestamp (as a whole number string)
function getTimeStamp() {
    return "" + round(getTime());
}

for (i = 0; i < imageList.length; i++) {
    // Select the original image and store its title
    selectWindow(imageList[i]);
    title = getTitle();

    // Find the position of ".lif - " in the title
    lifIndex = indexOf(title, ".lif - ");

    if (lifIndex > 0) {
        // Record the list of open windows before projection
        beforeList = getList("image.titles");

        // Perform Z projection on the current image stack
        run("Z Project...", "projection=[Max Intensity]");
        wait(500); // Wait for the projection window to appear

        // Get the list of open windows after projection
        afterList = getList("image.titles");

        // Find the new window by comparing afterList and beforeList
        projName = "";
        for (j = 0; j < afterList.length; j++) {
            found = false;
            for (k = 0; k < beforeList.length; k++) {
                if (afterList[j] == beforeList[k]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                projName = afterList[j];
                break;
            }
        }

        if (projName == "") {
            print("Error: Projection window not found for image: " + title);
            close(title);
            continue;
        }

        // Build a unique filename
        raw_image_name = substring(title, lifIndex + 7) + "_" + i + "_" + getTimeStamp() + ".tif";
        sanitized = replace(raw_image_name, " ", "_");
        sanitized = replace(sanitized, "/", "_");
        sanitized = replace(sanitized, "\\", "_");
        image_name = "batch" + batchNumber + "_" + sanitized;
        newFilename = outputFolder + image_name;
        print("Attempting to save: " + newFilename);

        // Save the projection window using the new filename
        selectWindow(projName);
        saveAs("Tiff", newFilename);
        print("Successfully saved: " + newFilename);

        // Close both the projection window and the original image
        close(projName);
        selectWindow(title);
        close();
    } else {
        close();
    }
}

// Close any remaining open images
run("Close All");
print("Batch " + batchNumber + " processing complete.");
