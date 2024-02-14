### Dependencies

# Basic utilities
import glob
import subprocess
import click

# Calculations
import numpy as np

# Image handling
from cv2 import imwrite

# Robust Photometric Stereo components
from rps import RPS # source repository script
import psutil # source repository script


### Shared functions

def list_images(img_prefix, img_dir="."):
    # Get the image file names of the files in the specified directory
    image_file_extensions = [
        "jpg",
        "jpeg",
        "JPG", 
        "JPEG", 
        "png", 
        "PNG", 
        "tif", 
        "tiff", 
        "TIF", 
        "TIFF"
    ]
    # Create a list of image files
    image_file_list = []
    for image_file_extension in image_file_extensions:
        image_file_list += (glob.glob(f"{img_dir}/{img_prefix}*.{image_file_extension}"))

    # Sort the file list to ensure the correct order of files
    image_file_list.sort()
    return(image_file_list)


def list_dirs(dir_prefix):
    stack_dir_list = glob.glob(f"{dir_prefix}*/")
    stack_dir_list.sort()
    return(stack_dir_list)

def image_file_list_checker(image_file_list):
    if len(image_file_list) == 0:
        print("No images to process in the specified directories")
        exit()


def write_image_file(rps, name=None):
    """
    Image exporter compatible with "RobustPhotometricStereo"
    """
    # Array to conventional image shape (height, width, 3 color channels)
    N = np.reshape(rps.N, (rps.height, rps.width, 3))
    # Swap RGB <-> BGR (BGR is the standard format of opencv)
    N[:, :, 0], N[:, :, 2] = N[:, :, 2], N[:, :, 0].copy()  
    # Rescale
    N = (N + 1.0) / 2.0  
    # Scale between 0 and 255
    N = (N - np.min(N)) / (np.max(N) - np.min(N)) * 255 
    # Convert to integers so most file formats can handle it
    N = np.uint8(N)
    # Write file
    imwrite(name, N)



### Shared Command line arguments

img_prefix_option = click.option(
    "-i",
    "--img_prefix",
    default="",
    help="Prefix of the desired image file names"
)

dir_prefix_option = click.option(
    "-d",
    "--dir_prefix",
    default="stack_",
    help="Prefix of the desired directory names"
)

mlic_option = click.option(
    "--mlic", 
    is_flag=True, 
    help="Shortcut for dealing with MLICs"
)


verbosity_option = click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity"
)

@click.group(
    help="A toolkit to process images recorded by focus stacking, panoramic imaging and multi light imaging", 
    epilog='Information on arguments and options of all commands are in the specific help messages: e.g. "focus-pano distribute --help"'
)
def main():
    pass



### Distribute images into folders of image stacks
@main.command(
    "distribute",
    short_help="Distribute stacks into directories"
)
@img_prefix_option
@dir_prefix_option
@mlic_option
@verbosity_option
@click.argument("stack_size", type=click.INT)


def distribute(img_prefix, dir_prefix, mlic, verbose, stack_size):
    """
    Distribute the specified images into directories (one focus stack or MLIC per directory)
    Specify the number of images in each focus stack/MLIC as an argument
    """

    if stack_size == None:
        print("Specify the size of the stacks with -s or --stack_size")
        exit()
    
    # Check MLIC mode
    if mlic:
        dir_name_base = "mlic_"
    else:
        dir_name_base = dir_prefix

    image_file_list = list_images(img_prefix)
 
    
    # Get the most important variables
    n_image_files = len(image_file_list)

    if verbose:
        print(f"Number of images found: {n_image_files}")
        print(f"Stack size: {stack_size}")
        print(f"Directory name prefix: {dir_name_base}")
    

    # Checking if the number of images fits the stack dimensions
    if n_image_files % stack_size != 0:
        print("Number of images does not match size of stack. Check for additional images like scales.")
        exit()

    # Calculate the number of stacks
    n_stacks = int(n_image_files / stack_size)


    # Create directories and remenber their names
    stack_dirs = []
    for i in range(1, n_stacks+1):
        numbered_dir = f"{dir_name_base}{str(i).zfill(3)}"
        subprocess.run(f"mkdir {numbered_dir}", shell=True)
        stack_dirs.append(numbered_dir)


    # Put the files into their corresponding directories
    image_file_counter = 0
    for i in range(0, n_stacks):
        image_list = " ".join(image_file_list[image_file_counter:(image_file_counter + stack_size)])
        subprocess.run(f"mv {image_list} {stack_dirs[i]}", shell=True)
        image_file_counter += stack_size



### Extract files from directories
@main.command(
    "dissolve",
    short_help='Dissolve directories (reverse of "distribute")'
)
@dir_prefix_option
@mlic_option


def dissolve(dir_prefix, mlic):
    """
    Extract all images from the specified directories into the working directory (reversing the "distribute" command)
    """

    # Check MLIC mode
    if mlic:
        dir_name_base = "mlic_"
    else:
        dir_name_base = dir_prefix

    # Create a list of directories 
    stack_dir_list = glob.glob(f"{dir_name_base}*/")

    for stack_dir in stack_dir_list:
        subprocess.run(f"mv {stack_dir}* .", shell=True)
        subprocess.run(f"rmdir {stack_dir}", shell=True)
 


### Merge focus stacks
@main.command(
    "focusmerge",
    short_help="Merge focus stacks",
    epilog="This command uses focus-stack (https://github.com/PetteriAimonen/focus-stack) for the alignment and either focus-stack or enfuse for merging"
)

@click.option(
    "-m",
    "--method",
    default="focus-stack",
    type=click.Choice(["focus-stack", "enfuse"]),
    help="Specifies the desired focus stack merging program"
)

@img_prefix_option
@dir_prefix_option
@mlic_option
@verbosity_option


def focusmerge(method, img_prefix, dir_prefix, mlic, verbose):
    """
    Focus merge all image files in the specified directories (matching a prefix, if specified)
    """
    # Create a list of directories 
    
    if verbose:
        print(f"Focus merging method: {method}")

    stack_dir_list = list_dirs(dir_prefix)

    for stack_dir in stack_dir_list:
        stack_dir_name = stack_dir[:-1]

        image_file_list = list_images(img_prefix, img_dir=stack_dir)
        
        if verbose:
            print(f"Processing directory {stack_dir} ({len(image_file_list)} images)")
        
        
        # Check if there are any image files to process
        image_file_list_checker(image_file_list)
        
        # Make the python list a long space delimited string 
        image_file_list_str = " ".join(image_file_list)
        
        # Do the actual focus merging for each directory
        if method == "focus-stack":
            subprocess.run(f"focus-stack {image_file_list_str} --no-contrast --no-whitebalance --output=merged_{stack_dir_name}_focus-stack.tif", shell=True)
        elif method == "enfuse":
            image_file_list_wo_dirs = " ".join([file_name.split("/")[1] for file_name in image_file_list])
            
            print(image_file_list_wo_dirs)

            subprocess.run(f"cd {stack_dir} && focus-stack --align-only --no-whitebalance {image_file_list_wo_dirs}", shell=True)
            subprocess.run(f"cd {stack_dir} && enfuse aligned_* --hard-mask --contrast-weight 1 --exposure-weight 0 --saturation-weight 0 --contrast-window 7 --contrast-edge-scale=0.3 --contrast-min-curvature=-0.5% -o merged_{stack_dir[:-1]}_enfuse.tif", shell=True)
            
            subprocess.run(f"cd {stack_dir} && mv *_enfuse.tif .. && rm aligned*", shell=True)



### Photometric Stereo
@main.command(
    "normalmapper",
    short_help="Perform photometric stereo",
    epilog="This command relies on RobustPhotometricStereo (https://github.com/yasumat/RobustPhotometricStereo) for creating normal maps. Note that the light point file has to be compliant with the format of Relight (https://github.com/cnr-isti-vclab/relight)"
)

@click.option(
    "-l",
    "--lp-file",
    default="*.lp",
    help="Light point file name/location"
)


@click.option(
    "--solver",
    default="l2",
    type=click.Choice(["l2", "rpca", "l1", "sbl"]),
    help='Specifies a photometric stereo method. "l2" is the original (Woodham 1980) and fastest method (not accounting for strongly reflective surfaces)'
)


@click.option(
    "-d",
    "--dir_prefix",
    default="mlic",
    help="Prefix of the desired directory names"
)

@verbosity_option
@img_prefix_option



def normalmapper(solver, lp_file, dir_prefix, img_prefix, verbose):
    """
    Perform photometric stereo, creating normal maps for all MLICs located in the specified directories
    """
    
    if verbose:
        print(f"Photometric stereo solver: {solver}")


    # Point to correct RPS solver
    if solver == "l2":
        ps_method = RPS.L2_SOLVER    # Least-squares # Woodham 1980
    elif solver == "rpca":
        ps_method = RPS.RPCA_SOLVER    # Robust PCA # Wu et al. 2010 robust photometric stereo via low-rank matrix completion and recover
    elif solver == "l1":
        ps_method = RPS.L1_SOLVER_MULTICORE    # L1 residual minimization # Takes too long
    elif solver == "sbl":
        ps_method = RPS.SBL_SOLVER_MULTICORE    # Sparse Bayesian Learning # Takes too long⎄
    
    # Get the number of images per MLIC from the lp file
    lp_file_name = glob.glob(lp_file)[0]
    lp_file = np.loadtxt(lp_file_name, skiprows=1, usecols=(1,2,3))
    mlic_size = lp_file.shape[0]


    # Create a list of the directories with the to process images
    mlic_dir_list = list_dirs(dir_prefix)


    for mlic_dir in mlic_dir_list:
        mlic_dir_name = mlic_dir[:-1]

        image_file_list = list_images(img_prefix, img_dir=mlic_dir_name)

        if verbose:
            print(f"Processing directory {mlic_dir_name} ({len(image_file_list)} images)")

        # Check if there are any image files to process
        image_file_list_checker(image_file_list)    

        # Check the correct number of images in the directory
        if len(image_file_list) % mlic_size != 0:
            print(f"Number of images in directory {stack_dir_name} in conflict with lp file")
            exit()
        
        # Instantiate Robust Photometric Stereo object
        rps = RPS()

        # Load light matrix
        rps.load_relight_lp(filename=lp_file_name)    

        # Load images
        # The directory name needs a slash, otherwise psutil's load_mages won't work
        rps.load_images(foldername=mlic_dir, ext=image_file_list[0].split(".")[-1])    # Load observations
        
        # Compute
        rps.solve(ps_method) 
        
        write_image_file(rps, f"normalmap_{solver}_{mlic_dir_name}.png")


### Grazing light images (TO DO)
# Extract images of one light direction from a MLIC and store it into folders that correspond to stacks

### HDR (TO DO)
# Just like the normalmapper but with enfuse with hdr settings


### Stitch panorama tiles (translation only)
@main.command(
    "stitch",
    short_help="Stitch translation panoramas",
    epilog="Like Hugin, this command uses PanoTools command line programs to create a panorama"
)

@verbosity_option
@click.option("-i", "--img_prefix", default="merged_", help="Prefix of the desired image file names")
@click.option("-c", "--crop-factor", type=click.FLOAT, default=1, help='Crop factor of the camera sensor relative to a "full frame" sensor (36 mm *24 mm)')
@click.option("--sensor_width", type=click.FLOAT, default=36, help="Width of the camera sensor in mm")
@click.option("-r", "--rotation", is_flag=True, help="Allows for rotation")
@click.option("-z", "--z-translation", is_flag=True, help="Allows for translation in z direction")
@click.option("-f", "--focal-length", type=click.INT, help="Focal length of the camera optics in mm. Do not use equivalent focal lengths")


def stitch(focal_length, sensor_width, crop_factor, img_prefix, rotation, z_translation, verbose):
    """
    Stitch a panoramic image from multiple images located in the working directory
    Only such transformations are allowed that correspond to mainly translation movements between the image/stack captures
    Use either the crop factor or the sensor width, not both unless you are aware of the consequences (crop factor no longer corresponding to "full frame")
    """

    # Calculating the field of view from the focal length and the sensor size
    hfov = 2 * np.arctan( (sensor_width / crop_factor) / focal_length) * (180/np.pi)
    
    if verbose:
        print(f"Focal length: {focal_length} mm")
        print(f"Sensor width: {sensor_width} mm")
        print(f"Crop factor: {crop_factor}")
        print(f"Horizontal field of view: {hfov} degrees")
        

    # Generate a Hugin project
    if verbose:
        print("Generating the Hugin project")
    subprocess.run(f"pto_gen -f {hfov} -o auto_trXY.pto {img_prefix}*", shell=True)
    

    # Generate control points
    if verbose:
        print("Searching control points")
    subprocess.run(f"cpfind  -o auto_trXY.pto auto_trXY.pto", shell=True)


    # Discard misleading control points
    if verbose:
        print("Discarding misleading control points")
    subprocess.run(f"cpclean -o auto_trXY.pto auto_trXY.pto", shell=True)
    
    # Set the allowed image transformations and write them to the project file
    # TrX and TrY are pure translations (without any rotation or else)
    # r is for allowing rotations
    # TrZ is for allowing translations in Z direction
    if rotation:
        if z_translation:
            subprocess.run(f"pto_var --opt=TrX,TrY,TrZ,r -o auto_trXYZr.pto auto_trXY.pto", shell=True)
        else:
            subprocess.run(f"pto_var --opt=TrX,TrY,r -o auto_trXYr.pto auto_trXY.pto", shell=True)
    else:
        if z_translation:
            subprocess.run(f"pto_var --opt=TrX,TrY,TrZ -o auto_trXYZ.pto auto_trXY.pto", shell=True)
        else:
            subprocess.run(f"pto_var --opt=TrX,TrY -o auto_trXY.pto auto_trXY.pto", shell=True)
    
    
    # Optimize the transformation parameters according to the specified parameters
    if verbose:
        print("Optimizing transformation parameters")
    subprocess.run(f"autooptimiser -n -o auto_trXY.pto auto_trXY.pto", shell=True)
    
    # Set the output properties
    # -p 0 stands for the rectilinear perspective
    # field of view and the canvas size are automatically calculated
    subprocess.run(f"pano_modify -p 0 --fov=AUTO --canvas=AUTO -o auto_trXY.pto auto_trXY.pto", shell=True)
    
    # Perform image transformations and export result
    if verbose:
        print("Performing image transformations, blending and exporting")
    subprocess.run(f"hugin_executor --stitching --prefix=hugin_trXY_pano_ auto_trXY.pto", shell=True)



if __name__=="__main__":
    main(prog_name="focus-pano")
