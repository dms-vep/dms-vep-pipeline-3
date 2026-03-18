export async function parseVegaSpecFromHTML(htmlContent) {
    // Use the browser's DOMParser to parse the HTML content
    const parser = new DOMParser();
    const doc = parser.parseFromString(htmlContent, "text/html");
  
    // Locate the script tag containing the Vega spec
    const scripts = doc.querySelectorAll("script");
    let vegaSpecString = "";
  
    scripts.forEach((script) => {
      if (script.textContent.includes("var spec =")) {
        vegaSpecString = script.textContent;
      }
    });
  
    if (!vegaSpecString) {
      throw new Error("Vega spec not found in the HTML document.");
    }
  
    // Locate the JSON part after "var spec ="
    const specStart = vegaSpecString.indexOf("var spec =") + "var spec =".length;
    const specEnd = vegaSpecString.indexOf("};", specStart) + 1; // +1 to include the closing brace
    const jsonSubstring = vegaSpecString.substring(specStart, specEnd);
  
    // Trim any leading or trailing non-JSON characters (like semicolons or whitespace)
    const trimmedJsonString = jsonSubstring.trim();
  
    // Parse the JSON string into an object
    let vegaSpecObject;
    try {
      vegaSpecObject = JSON.parse(trimmedJsonString);
    } catch (error) {
      console.log(trimmedJsonString);
      throw new Error("Error parsing Vega spec JSON: " + error.message);
    }
  
    return vegaSpecObject;
  }