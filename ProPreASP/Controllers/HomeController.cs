using Microsoft.AspNetCore.Mvc;
using ProPreASP.Models;
using System.Diagnostics;
using System.Net;

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


        [HttpPost]
        [ValidateAntiForgeryToken]
        public async Task<IActionResult> dna(string dnaSeq)
        {
            // ATGGGTCCTTCAGTAGTTCCTATTAACCCA
            var url = "http://127.0.0.1:56000/GetAmino/" + dnaSeq;
            var userSessionId = "a" + HttpContext.Session.Id;
            string webRootPath = _webHostEnvironment.WebRootPath;
            /*
             GACCCGGCAGGCCTTGCGCGGGCAACATGGCGGCGCCCGGCGAGCGGGGCCGCTTCCACGGCGGGAACCTCTTCTTCCTGCCGGGGGGCGCGCGCTCCGAGATGATGGACGACCTGGCGACCGACGCGCGGGGCCGGGGCGCGGGGCGGAGAGACGCGGCCGCCTCGGCCTCGACGCCAGCCCAGGCGCCGACCTCCGATTCTCCTGTCGCCGAGGACGCCTCCCGGAGGCGGCCGTGCCGGGCCTGCGTCGACTTCAAGACGTGGATGCGGACGCAGCAGAAGGTGCAGTTCCCTGCCCGATTTCTCCCAGCCCCGCGCAGCCCCTGTCCCCGCCCCCGCCCAGGTACCCCGGCAGAGCTTCCCAGGGTTGCCTGTCCCTGAACCTTGCCCCCCGGGTAGGCCCGGCCTTACAGCCTTCATCCGCGCGTGGGTTGGATCGTCTGCAGGACTTTGGCCGGAGTCCAGTGGGCCACCGGCTGGGCCGTACAGTGGGGAGCTTTGGGCGCCTTTGTTCGGAGAATGAACTCACTCTCGGTCGGCCTGCTTCCGCAGCGGGAC
             */
            await Task.Run(() =>
             DownloadFile(url, webRootPath + "\\pdbFiles", userSessionId + ".pdb"));

            var filePath = userSessionId + ".pdb";
            return RedirectToAction("dna3d", new {filePath});
        }

        public IActionResult dna3d(string filePath)
        {
            ViewData["filePath"] = filePath;
            return View();
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
            using (var result = await client.GetAsync(url))
                return result.IsSuccessStatusCode ? await result.Content.ReadAsByteArrayAsync() : null;
        }


        public IActionResult protien()
        {
            return View();
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