# Provide methods for accessing imgCIF data

#==

# Background

imgCIF is a text format for describing raw data images collected on
crystallographic instruments. As such it is not suited for actually
containing the raw data, but instead provides tags pointing to the
data storage location.

This library provides methods for bringing that image data into
Julia, including handling the variety of formats provided, and for
general manipulation of the information found in the imgCIF.

==#
"""
Handle items in an imgCIF file. See method `imgload`.
"""
module ImgCIFHandler

using CrystalInfoFramework
using DataFrames
using InvertedIndices
import Tar
using CodecBzip2
using CodecZlib
using CodecInflate64    #for Windows Zip files
using Downloads: Curl, Downloader, download, Response, RequestError, request
using ZipArchives

# See format-specific includes for more using statements

using TranscodingStreams
using URIs
using SimpleBufferStream
using LinearAlgebra
using ImageBinarization
using ImageFiltering
using Statistics
using Rotations

export create_archives  #register an archive
export imgload         #Load raw data
export peek_image      #Find first image in archive
export ping_archive    #Check that URL exists
export make_absolute_uri #Use Cif block contents to make absolute URI
export get_detector_axis_settings #Get axis settings for particular frame
export get_beam_centre
export get_gonio_axes
export get_pixel_coordinates
export get_detector_distance   #Get distance of flat detector from crystal
export get_surface_axes #get the axes used to locate pixels on the detector
export get_id_sequence #get a list of sequential binary ids from the same scan
export get_axis_vector #get the vector for an axis_id
export find_peaks #Find some peaks in an image
export bin_id_from_scan_frame #Convert scan/frame_no into binary id
export scan_frame_from_bin_id #Convert bin_id into scan/frame
export peak_to_frames #Calculate all peak appearances
export get_dependency_chain #Get the dependent axes of an axis

export Peak
export intensity, coords, frame, scan, dist #Working with peaks

export ImageArchive
export has_local_version #for ImageArchive types
export get_constant_part #for RsyncArchive work

include("imgcif_base.jl")
include("hdf_image.jl")
include("cbf_image.jl")
include("adsc_image.jl")
include("kcd_image.jl")
include("imgcif.jl")
include("recip.jl")
include("zipfilemanager.jl")


end
