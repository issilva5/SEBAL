package org.fogbowcloud.sebal.parsers.plugins.ftp;

import static org.mockito.Mockito.doReturn;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.spy;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Properties;

import org.apache.commons.io.FileUtils;
import org.apache.http.HttpException;
import org.fogbowcloud.sebal.parsers.plugins.StationOperatorConstants;
import org.json.JSONArray;
import org.json.JSONObject;
import org.junit.Assert;
import org.junit.Test;

public class TestFTPStationOperator {
	
	@Test
	public void testFindNearestStationCorrectCalculation() throws URISyntaxException, HttpException, IOException, ParseException {
		
		// set up
		Properties properties = mock(Properties.class);
		String year = "2002";
		int numberOfDays = 0;
		double lat = -3.40;
		double lon = -45.20;
		
		String stringDate = "26-01-2002";		
		SimpleDateFormat f = new SimpleDateFormat("dd-MM-yyyy");
		Date d = f.parse(stringDate);
		long milliseconds = d.getTime();
		Date date = new Date(milliseconds);
		
		JSONObject stationOne = new JSONObject();
		stationOne.put("id", "82294");
		stationOne.put("lon", "-40.13333333");
		stationOne.put("altitude", "16.5");
		stationOne.put("name", "ACARAU - CE");
		stationOne.put("lat", "-2.88333333");
		
		JSONObject stationTwo = new JSONObject();
		stationTwo.put("id", "83096");
		stationTwo.put("lon", "-37.05");
		stationTwo.put("altitude", "4.72");
		stationTwo.put("name", "ARACAJU - SE");
		stationTwo.put("lat", "-10.95");
		
		JSONArray stations = new JSONArray();
		stations.put(stationOne);
		stations.put(stationTwo);
		
		List<JSONObject> expectedStation = new ArrayList<JSONObject>();
		expectedStation.add(stationOne);
		
		FTPStationOperator stationOperator = spy(new FTPStationOperator(properties));
		doReturn(stations).when(stationOperator).getStations(year);
		
		// exercise
		List<JSONObject> chosenStation = stationOperator.findNearestStation(date, lat, lon, numberOfDays);
		
		// expect
		Assert.assertEquals(expectedStation, chosenStation);
	}
	
	@Test
	public void testReadStation() throws Exception {
		// set up
		String stringDate = "26-01-2002";		
		SimpleDateFormat f = new SimpleDateFormat("dd-MM-yyyy");
		Date d = f.parse(stringDate);
		long milliseconds = d.getTime();
		Date date = new Date(milliseconds);		
		
		// Establishing local fixed attributes
		String year = "2002";
		String compressedUnformattedStationFilePath = "/tmp/2002/827910-99999-2002.gz";
		String stationFileUrl = "fake-station-file-url";
		File compressedUnformattedStationFile = new File(compressedUnformattedStationFilePath);
		Properties properties = mock(Properties.class);

		// Creating temporary year directory
		String baseUnformattedLocalStationFilePath = "/tmp/2002";
		File baseUnformattedLocalStationDir = new File(baseUnformattedLocalStationFilePath);
		baseUnformattedLocalStationDir.mkdirs();
		
		// Copying station file to be used from resource to /tmp
		String originalCompressedFilePath = "src/test/resource/827910-99999-2002.gz";
		File originalCompressedFile = new File(originalCompressedFilePath);
		FileUtils.copyFile(originalCompressedFile, compressedUnformattedStationFile);
		
		FTPStationOperator stationOperator = spy(new FTPStationOperator(properties));
		doReturn(baseUnformattedLocalStationFilePath).when(stationOperator).getBaseUnformattedLocalStationFilePath(year);
		doReturn(compressedUnformattedStationFile).when(stationOperator).getUnformattedStationFile("82791", year);
		doReturn(stationFileUrl).when(stationOperator).getStationFileUrl("82791", year);
		doReturn(true).when(stationOperator).downloadUnformattedStationFile(compressedUnformattedStationFile, stationFileUrl);
		
		// exercise
		JSONArray stationData = stationOperator.readStation("82791",
				StationOperatorConstants.DATE_FORMAT.format(date),
				StationOperatorConstants.DATE_FORMAT.format(date));
		
		// expect
		Assert.assertNotNull(stationData);
	}
	
	@Test
	public void testGetStations() throws IOException {
		// set up
		Properties properties = new Properties();
		properties.put("station_ftp_server_url", "http://www2.lsd.ufcg.edu.br/~esdras"); // TODO: transfer station csv file from resource to this path
		
		String originalStationCSVFilePath = "src/test/resource/2002-stations.csv";
		File originalStationCSVFile = new File(originalStationCSVFilePath);
		
		String destinyStationCSVFilePath = "/tmp/2002/2002-stations.csv";
		File destinyStationCSVFile = new File(destinyStationCSVFilePath);
		
		FileUtils.copyFile(originalStationCSVFile, destinyStationCSVFile);
		
		String year = "2002";
		String localStationsCSVFilePath = "/tmp/2002";
		File localStationsCSVDir = new File(localStationsCSVFilePath);
		localStationsCSVDir.mkdirs();
		
		FTPStationOperator stationOperator = spy(new FTPStationOperator(properties));
		doReturn(localStationsCSVFilePath).when(stationOperator).getStationCSVDirPath(year);
		
		// exercise
		JSONArray stations = stationOperator.getStations(year);
		FileUtils.deleteDirectory(localStationsCSVDir);
		
		// expect
		Assert.assertNotNull(stations);
	}
}
