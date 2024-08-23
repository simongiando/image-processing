//** Only change this
var roi = ee.FeatureCollection('---------------')
var region_name = '------'

//** Downloading collection
var year_one = 2014
var year_ten = 2024

var start_date = 0
var end_date = 365

//** Downloading one year
var year_start = 2023
var year_end = 2023

var month_start = 1
var month_end = 12

var day_start = 1
var day_end = 31

//Source: https://code.earthengine.google.com/c0074d0939db1de9d0cd04673f577ccb?accept_repo=users%2Femaprlab%2Fpublic
// Import, filter, mask, and join all the landsats imagery. Simplified for sanity's sake.
function load_landsat_imagery (start_year, end_year,start_day, end_day, filter_geometry, input_mask){
    //Define the masking fuction, masking is more efficent than clipping
    function apply_mask(image){
      return ee.Image(image).updateMask(input_mask).toUint16();
    }

    //Load the Landsat Imagery
    var landsat_4 = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2') 
        .filterBounds(filter_geometry)
        .filterDate(ee.Date.fromYMD(start_year,1,1), ee.Date.fromYMD(end_year,12,31)) 
        .filter(ee.Filter.calendarRange(start_day, end_day))
        .select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','QA_PIXEL'],['B1','B2','B3','B4','B5','B7','pixel_qa']) 
        .map(apply_mask);
    var landsat_5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2') 
        .filterBounds(filter_geometry)
        .filterDate(ee.Date.fromYMD(start_year,1,1), ee.Date.fromYMD(end_year,12,31)) 
        .filter(ee.Filter.calendarRange(start_day, end_day))
        .select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','QA_PIXEL'],['B1','B2','B3','B4','B5','B7','pixel_qa']) 
        .map(apply_mask);
    var landsat_7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2') 
        .filterBounds(filter_geometry)
        .filterDate(ee.Date.fromYMD(start_year,1,1), ee.Date.fromYMD(end_year,12,31)) 
        .filter(ee.Filter.calendarRange(start_day, end_day))
        .select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','QA_PIXEL'],['B1','B2','B3','B4','B5','B7','pixel_qa']) 
        .map(apply_mask);
    var landsat_8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') 
        .filterBounds(filter_geometry)
        .filterDate(ee.Date.fromYMD(start_year,1,1), ee.Date.fromYMD(end_year,12,31)) 
        .filter(ee.Filter.calendarRange(start_day, end_day))
        .select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL'],['B1','B2','B3','B4','B5','B7','pixel_qa']) 
        .map(oli_to_etm)
        .map(apply_mask);
    var landsat_9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2') 
        .filterBounds(filter_geometry)
        .filterDate(ee.Date.fromYMD(start_year,1,1), ee.Date.fromYMD(end_year,12,31)) 
        .filter(ee.Filter.calendarRange(start_day, end_day))
        .select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL'],['B1','B2','B3','B4','B5','B7','pixel_qa']) 
        .map(oli_to_etm)
        .map(apply_mask);
    
    var merged = ee.ImageCollection(landsat_9.merge(landsat_8).merge(landsat_7).merge(landsat_4).merge(landsat_5));

    return merged;
}

function apply_pixel_qa (img){
    var qa = img.select('pixel_qa');
    var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 5).eq(0));
    return img.updateMask(mask);
}

//SOURCE: Emapr lab public landtrendr module
//https://code.earthengine.google.com/c0074d0939db1de9d0cd04673f577ccb?accept_repo=users%2Femaprlab%2Fpublic
//------ L8 to L7 HARMONIZATION FUNCTION -----
// slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhgang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, Remote Sensing of Environment, 185, 57-70.(http://dx.doi.org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients
function oli_to_etm (oli){
    //Convert OLI reflectance to ETM using the Roy et al. 2018 RMA coefficents
    var band = oli.select(['B1','B2','B3','B4','B5','B7']);
    //RMA - create an image of slopes per band for L8 TO L7 regression line - David Roy

    var slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949]);
    var itcp = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029]);
    var y = band.resample('bicubic').subtract(itcp.multiply(10000)).divide(slopes)
        .set('system:time_start', oli.get('system:time_start'));
    
    return y.addBands(oli.select('pixel_qa')).toInt16();
}

//SOURCE: Emapr lab public landtrendr module
//https://code.earthengine.google.com/c0074d0939db1de9d0cd04673f577ccb?accept_repo=users%2Femaprlab%2Fpublic
// mediod function for annual compositing
function medoid_mosaic (inCollection){
    //Create a medoid composite from the collection of input Landsat images'''
    //fill in missing years with the dummy collection
    var dummyCollection = ee.ImageCollection([
        ee.Image([0,0,0,0,0,0]).mask(ee.Image(0)).rename(['B1','B2','B3','B4','B5','B7'])]);

    var finalCollection = inCollection.merge(dummyCollection);
    
    //calculate median across images in collection per band
    var median = finalCollection.median();                                 

    //calculate the different between the median and the observation per image per band
    function inner_map (img){
        var diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));                                  
        return diff.reduce('sum').addBands(img);
    }
    var difFromMedian = finalCollection.map(inner_map);
    
    return ee.ImageCollection(difFromMedian).reduce(ee.Reducer.min(7)).select([1,2,3,4,5,6], ['B1','B2','B3','B4','B5','B7']);
}

function create_annual_composites(input_collection, start_year, end_year){
    var composited_images = [];
    for (var year = start_year; year <= end_year; year++) {
        var start_date = ee.Date.fromYMD(year, 1, 1);
        var end_date = ee.Date.fromYMD(year, 12, 31);
        var composite = medoid_mosaic(input_collection.filterDate(start_date, end_date)) 
            .set('system:time_start', ee.Date.fromYMD(year, 8, 1).millis()) 
            .toFloat();
         
        composited_images = composited_images.concat(composite);
    }
    return ee.ImageCollection.fromImages(composited_images);
}
var calculate_nbr = function(image){

    var nbr = image.normalizedDifference(["B4", "B7"]).rename(["NBR"]).multiply(-1);

    return nbr.addBands(image)
        .set("system:time_start", image.get("system:time_start"))
        .float();
    
};

//########################################################################################################
//##### UNPACKING LT-GEE OUTPUT STRUCTURE FUNCTIONS ##### 
//########################################################################################################
//Source: https://code.earthengine.google.com/c0074d0939db1de9d0cd04673f577ccb?accept_repo=users%2Femaprlab%2Fpublic
// ----- FUNCTION TO EXTRACT VERTICES FROM LT RESULTS AND STACK BANDS -----
var flatten_ltr_fits = function(ltr_output, start_year, end_year){

    start_year = ee.Number(start_year);
    end_year = ee.Number(end_year);

    var years = ee.List.sequence(start_year, end_year)  ;  
    
    var getYearStr = function(year){
          return(ee.String('yr_').cat(ee.Algorithms.String(year).slice(0,4)))};

    var yearsStr = years.map(getYearStr);
                        
    var fitted_values = ltr_output.select(['B1_fit','B2_fit','B3_fit','B4_fit','B5_fit','B7_fit'])
                     .arrayFlatten([yearsStr])
                     .toInt16();
                        
    var b1_fits = ltr_output.select(['B1_fit']).arrayFlatten([yearsStr]).toInt16();
    var b2_fits = ltr_output.select(['B2_fit']).arrayFlatten([yearsStr]).toInt16();
    var b3_fits = ltr_output.select(['B3_fit']).arrayFlatten([yearsStr]).toInt16();
    var b4_fits = ltr_output.select(['B4_fit']).arrayFlatten([yearsStr]).toInt16();
    var b5_fits = ltr_output.select(['B5_fit']).arrayFlatten([yearsStr]).toInt16();
    var b7_fits = ltr_output.select(['B7_fit']).arrayFlatten([yearsStr]).toInt16()  ;               
                        
    var reshuffle = function(year){
        var image = b1_fits.select([year]).rename(['LS_B'])
            .addBands(b2_fits.select([year]).rename(['LS_G']))
            .addBands(b3_fits.select([year]).rename(['LS_R']))
            .addBands(b4_fits.select([year]).rename(['LS_NIR']))
            .addBands(b5_fits.select([year]).rename(['LS_SWIR1']))
            .addBands(b7_fits.select([year]).rename(['LS_SWIR2']));
        return image.toInt16();
    };

    var all_images = ee.ImageCollection.fromImages(yearsStr.map(reshuffle));

    return all_images;
};

// Runs landtrendr over the annual composites using NBR
//SOURCE: https://emapr.github.io/LT-GEE/running-lt-gee.html
var interpolate_annual_composites = function(image_collection, start_year, end_year){

    image_collection = image_collection.map(calculate_nbr);

    var ltr_output = ee.Algorithms.TemporalSegmentation.LandTrendr({
        timeSeries:image_collection,
        maxSegments:9,
        spikeThreshold:0.9,
        vertexCountOvershoot:3,
        preventOneYearRecovery:true,
        recoveryThreshold:0.25,
        pvalThreshold:0.05,
        bestModelProportion:0.75,
        minObservationsNeeded:6
    });

    var flattend_images = flatten_ltr_fits(ltr_output, start_year, end_year);

    return flattend_images;
};

//Assigns a date to the new annual composites, corresponding to August 1st of each year (roughly the middle)
var set_date = function(col_bands,start_year, end_year){
  var imgs = [];                                                                    // create empty array to fill
  for (var i = start_year; i <= end_year; i++) {
    var year_string = i - start_year ;
    year_string = year_string.toString().concat('_.*');
    var year_img = col_bands.select(year_string);
    year_img = year_img.set('system:time_start', ee.Date.fromYMD(i, 8, 1).millis());
    imgs = imgs.concat(year_img.rename(['B1','B2','B3','B4','B5','B7']));
    }
  return ee.ImageCollection.fromImages(imgs);
};

var generate_historical_composites = function(tile, start_year, end_year, start_day, end_day){
    //create a buffer for invariant pixel selection (NOT USED IN THIS SCRIPT)
    var geometry_buffered = roi
    
    // Get a mask of the study area and the buffered area (for later clipping)
    var study_area_mask_buffer = ee.Image.constant(0).byte().paint(geometry_buffered, 1).byte();

    //Load in the cross-harmonized, temporally fileted Landsat time-series
    var raw_imagery = load_landsat_imagery(start_year, end_year,start_day, end_day, geometry_buffered, study_area_mask_buffer);
    
    //Apply qa to remove clouds
    raw_imagery = raw_imagery.map(apply_pixel_qa).select(['B1','B2','B3','B4','B5','B7']);
    
    //Produce the annual composites using the medoid methodology
    var medoids = create_annual_composites(raw_imagery, start_year, end_year);
    
    //Interpolate annual composites using Landtrendr
    var interpolated_medoids = interpolate_annual_composites(medoids, start_year, end_year);
    
    //Seperate images for each year, and assign them a date
    interpolated_medoids = set_date(interpolated_medoids.toBands(),start_year, end_year);

    return interpolated_medoids;
};

//** Generate composite
var composites = generate_historical_composites(roi,year_one,year_ten, start_date, end_date);

//** Generate composite only one year
var one_year = composites.filterDate(ee.Date.fromYMD(year_start, month_start, day_start), ee.Date.fromYMD(year_end, month_end, day_end)).median();

Map.addLayer(one_year.select(['B3','B2','B1']), {min:5000, max:20000});

Export.image.toAsset({
  image: one_year, 
  description: region_name+'_Landsat_'+year_start,
  region: roi,
  scale: 30,
  maxPixels: 1e13,
})
