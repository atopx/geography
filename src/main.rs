use std::f64::consts::PI;

// Pi*3000/180
const XPI: f64 = 52.35987755982988;
// 长半轴
const AXIS: f64 = 6378245.0;
// 偏心率平方
const OFFSET: f64 = 0.00669342162296594323;
// 百度坐标系经度偏移
const BD_LON_OFFSET: f64 = 0.0065;
// 百度坐标系纬度偏移
const BD_LAT_OFFSET: f64 = 0.0060;
const BD_Z_OFFSET: f64 = 0.00002;
const BD_THETA_OFFSET: f64 = 0.000003;
const ACCURACY_THRESHOLD: f64 = 0.0000000001;
const WGS_OFFSET: f64 = 0.01;

fn in_china(lon: f64, lat: f64) -> bool {
    (73.66 < lon && lon < 135.05) && (3.86 < lat && lat < 53.55)
}


fn convert(lon: f64, lat: f64) -> (f64, f64) {
    let lonlat = lon * lat;
    let abs_x = lon.abs().sqrt();
    let lon_pi = lon * PI;
    let lat_pi = lat * PI;
    let d = (lon_pi * 6.0).sin() * 20.0 + (lon_pi * 2.0).sin() * 20.0;
    let mut x = d;
    let mut y = d;
    x += 20.0 * lat_pi.sin() + (lat_pi / 3.0).sin() * 40.0;
    x += 160.0 * (lat_pi / 12.0).sin() + 320.0 * (lat_pi / 30.0).sin();

    y += 20.0 * lon_pi.sin() + (lon_pi / 3.0).sin() * 40.0;
    y += 150.0 * (lon_pi / 12.0).sin() + 300.0 * (lon_pi / 30.0).sin();
    let threshold: f64 = 2.0 / 3.0;
    x *= threshold;
    y *= threshold;
    x += 2.0 * lon + 3.0 * lat + 0.2 * lat * lat + 0.1 * lonlat + 0.2 * abs_x - 100.0;
    y += lon + 2.0 * lat + 0.1 * lon * lon + 0.1 * lonlat + 0.1 * abs_x + 300.0;
    return (x, y);
}

fn delta(lon: f64, lat: f64) -> (f64, f64) {
    let (de_lat, de_lon) = convert(lon - 105.0, lat - 35.0);
    let rad_lat = lat / 180.0 * PI;
    let magic = rad_lat.sin();
    let magic = 1.0 - OFFSET * magic * magic;
    let sqrtmagic = magic.sqrt();
    let de_lat = (de_lat * 180.0) / ((AXIS * (1.0 - OFFSET)) / (magic * sqrtmagic) * PI);
    let de_lon = (de_lon * 180.0) / (AXIS / sqrtmagic * rad_lat.cos() * PI);
    return (lon + de_lon, lat + de_lat);
}


fn wgs84_to_gcj02(lon: f64, lat: f64) -> (f64, f64) {
    if !in_china(lon, lat) {
        return (lon, lat);
    }
    return delta(lon, lat);
}

fn gcj02_to_bd09(lon: f64, lat: f64) -> (f64, f64) {
    let z = (lon * lon + lat * lat).sqrt() + BD_Z_OFFSET * (lat * XPI).sin();
    let theta = lat.atan2(lon) + BD_THETA_OFFSET * (lon * XPI).cos();
    return (z * theta.cos() + BD_LON_OFFSET, z * theta.sin() + BD_LAT_OFFSET);
}

fn bd09_to_gcj02(lon: f64, lat: f64) -> (f64, f64) {
    let x = lon - BD_LON_OFFSET;
    let y = lat - BD_LAT_OFFSET;
    let z = (x * x + y * y).sqrt() - BD_Z_OFFSET * (y * XPI).sin();
    let theta = y.atan2(x) - BD_THETA_OFFSET * (x * XPI).cos();
    return (z * theta.cos(), z * theta.sin());
}

fn gcj02_to_wgs84(lon: f64, lat: f64) -> (f64, f64) {
    let mut mlon = lon - WGS_OFFSET;
    let mut mlat = lat - WGS_OFFSET;
    let mut plon = lon + WGS_OFFSET;
    let mut plat = lat + WGS_OFFSET;
    let mut dlon: f64;
    let mut dlat: f64;
    let mut wgs_lon: f64;
    let mut wgs_lat: f64;
    loop {
        wgs_lat = (mlat + plat) / 2.0;
        wgs_lon = (mlon + plon) / 2.0;
        let (tmp_lon, tmp_lat) = delta(wgs_lon, wgs_lat);
        dlon = tmp_lon - wgs_lon;
        dlat = tmp_lat - wgs_lat;
        if dlat.abs() < ACCURACY_THRESHOLD && dlon.abs() < ACCURACY_THRESHOLD {
            break;
        }
        if dlat > 0.0 {
            plat = wgs_lat;
        } else {
            mlat = wgs_lat;
        }
        if dlon > 0.0 {
            plon = wgs_lon;
        } else {
            mlon = wgs_lon;
        }
    }
    return (wgs_lon, wgs_lat);
}


#[test]
fn wgs84_to_gcj02_test() {
    assert_eq!(
        wgs84_to_gcj02(116.65383912049425, 39.992875254426366),
        (116.6597270000643, 39.99397599992571)
    )
}

#[test]
fn gcj02_to_bd09_test() {
    assert_eq!(
        gcj02_to_bd09(113.56733, 34.81754),
        (113.5739274927449, 34.82327621798387)
    )
}

#[test]
fn bd09_to_gcj02_test() {
    assert_eq!(
        bd09_to_gcj02(113.5739274927449, 34.82327621798387),
        (113.56732983450783, 34.81754111708383)
    )
}

#[test]
fn gcj02_to_wgs84_test() {
    assert_eq!(
        gcj02_to_wgs84(113.5739274927449, 34.82327621798387),
        (113.56784261660992, 34.82437523040346)
    )
}

fn main() {
    let (x, y) = bd09_to_gcj02(113.5739274927449, 34.82327621798387);
    println!("{}, {}", x, y);
}
