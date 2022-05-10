using Microsoft.AspNetCore.Mvc;
using ProPreASP.Models;
using System.Diagnostics;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

namespace ProPreASP.Controllers
{
    public class HomeController : Controller
    {
        private readonly ILogger<HomeController> _logger;
        private readonly IWebHostEnvironment _webHostEnvironment;

        public HomeController(ILogger<HomeController> logger, IWebHostEnvironment webHostEnvironment)
        {
            _logger = logger;
            _webHostEnvironment = webHostEnvironment;
        }

        public IActionResult Index()
        {
            return View();
        }

        public IActionResult dna()
        {
            return View();
        }

        public IActionResult Embedded3D(string filePath)
        {
            ViewData["filePath"] = filePath;
            return View();
        }

        public async Task<JObject> GetJsonAlign(string dnaSeq, string diseaseName)
        {
            using (HttpClient hc = new HttpClient())
            {
                hc.Timeout = TimeSpan.FromSeconds(900);
                var json = (await hc.GetStringAsync($"http://localhost:57000/GetAlignScore/{dnaSeq}/{diseaseName}"))
                    .Replace("\u0022", "\"").Replace("\n", "");
                return JObject.Parse(json);
            }
        }

        

        [HttpPost]
        [ValidateAntiForgeryToken]
        public async Task<IActionResult> dna(string dnaSeq, string diseaseName, IFormFile file)
        {
            // ATGGGTCCTTCAGTAGTTCCTATTAACCCA
            var url = "http://127.0.0.1:57000/GetAminoSeq/" + dnaSeq;
            string webRootPath = _webHostEnvironment.WebRootPath;
            if (file != null)
            {
                using (var ms = new MemoryStream())
                {
                    file.CopyTo(ms);
                    var fileBytes = ms.ToArray();
                    dnaSeq = "";
                    if (!Directory.Exists(webRootPath + "\\pdbFiles")) Directory.CreateDirectory(webRootPath + "\\pdbFiles");
                    await System.IO.File.WriteAllBytesAsync($"{webRootPath + "\\pdbFiles"}/{file.FileName}", fileBytes);
                    foreach (var item in await System.IO.File.ReadAllLinesAsync($"{webRootPath + "\\pdbFiles"}/{file.FileName}"))
                    {
                        if (!item.StartsWith('>'))
                        {
                            dnaSeq += item;
                        }
                    }
                    url = "http://127.0.0.1:57000/GetAminoSeq/" + dnaSeq;
                }
            }
            var userSessionId = "a" + HttpContext.Session.Id;
            await Task.Run(() =>
             DownloadFile(url, webRootPath + "\\pdbFiles", userSessionId + ".pdb"));

            var filePath = userSessionId + ".pdb";
            return RedirectToAction("dna3d", new {filePath, dnaSeq, diseaseName});
        }

        public async Task<IActionResult> dna3d(string filePath, string dnaSeq, string diseaseName)
        {
            if (diseaseName == "rheumatoid arthritis cardiac")
            {
                var data = await GetJsonAlign(dnaSeq, "Homo-il34-v1.fasta");
                ViewData["filePath"] = filePath;
                ViewData["Score"] = Convert.ToDouble(data["Score"]);
                ViewData["Alignment1"] = data["Alignment"][0];
                ViewData["Alignment2"] = data["Alignment"][1];
                return View();
            }
            else if (diseaseName == "rheumatoid arthritis in bones")
            {
                var data = await GetJsonAlign(dnaSeq, "Homo Sapiens RD44 1.fasta");
                ViewData["filePath"] = filePath;
                ViewData["Score"] = Convert.ToDouble(data["Score"]);
                ViewData["Alignment1"] = data["Alignment"][0];
                ViewData["Alignment2"] = data["Alignment"][1];
                return View();
            }
            return BadRequest();
            //return Redirect("~/index.html");
        }

        static async Task DownloadFile(string url, string pathToSave, string fileName)
        {
            var content = await GetUrlContent(url);
            if (content != null)
            {
                if (!Directory.Exists(pathToSave)) Directory.CreateDirectory(pathToSave);
                await System.IO.File.WriteAllBytesAsync($"{pathToSave}/{fileName}", content);
            }
        }

        static async Task<byte[]?> GetUrlContent(string url)
        {
            using (var client = new HttpClient())
            {
                client.Timeout = TimeSpan.FromSeconds(900);
                using (var result = await client.GetAsync(url))
                    return result.IsSuccessStatusCode ? await result.Content.ReadAsByteArrayAsync() : null;
            }
        }


        public IActionResult protien()
        {
            return View();
        }


        [HttpPost]
        [ValidateAntiForgeryToken]
        public async Task<IActionResult> Protein(IFormFile file, string diseaseName)
        {
            try
            {
                if (file != null && file.Length > 0)
                {
                    using (var client = new HttpClient())
                    {
                        try
                        {
                            client.BaseAddress = new Uri("http://127.0.0.1:57000/PDB2SEQ");

                            byte[] data;
                            using (var br = new BinaryReader(file.OpenReadStream()))
                                data = br.ReadBytes((int)file.OpenReadStream().Length);

                            ByteArrayContent bytes = new ByteArrayContent(data);


                            MultipartFormDataContent multiContent = new MultipartFormDataContent();

                            multiContent.Add(bytes, "file", file.FileName);
                            using (var result = await client.PostAsync(client.BaseAddress, multiContent))
                            {
                                var newResult = System.Text.Encoding.UTF8.GetString(await result.Content.ReadAsByteArrayAsync());
                                string dnaSeq = newResult;
                                var url = "http://127.0.0.1:57000/GetAminoSeq/" + newResult;
                                var userSessionId = "a" + HttpContext.Session.Id;
                                string webRootPath = _webHostEnvironment.WebRootPath;
                                await Task.Run(() =>
                                 DownloadFile(url, webRootPath + "\\pdbFiles", userSessionId + ".pdb"));

                                var filePath = userSessionId + ".pdb";
                                return RedirectToAction("dna3d", new { filePath, dnaSeq, diseaseName });
                            }
                        }
                        catch (Exception)
                        {
                            return StatusCode(500); // 500 is generic server error
                        }
                    }
                }

                return StatusCode(400); // 400 is bad request

            }
            catch (Exception)
            {
                return StatusCode(500); // 500 is generic server error
            }
        }


        public IActionResult Privacy()
        {
            return View();
        }

        [ResponseCache(Duration = 0, Location = ResponseCacheLocation.None, NoStore = true)]
        public IActionResult Error()
        {
            return View(new ErrorViewModel { RequestId = Activity.Current?.Id ?? HttpContext.TraceIdentifier });
        }
    }
}