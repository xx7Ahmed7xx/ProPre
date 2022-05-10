using System.Text;

namespace ProPreASP.Controllers
{
    public static class HelperClass
    {
        /// <summary>
        /// Reads a form file from request body and returns it's text content.
        /// </summary>
        /// <param name="file">The file object from the form.</param>
        /// <returns></returns>
        public async static Task<string> ReadFastaAsync(this IFormFile file)
        {
            var result = new List<string>();
            using (var reader = new StreamReader(file.OpenReadStream()))
            {
                while (reader.Peek() >= 0)
                    result.Append(await reader.ReadLineAsync());
            }
            var final = "";
            foreach (var item in result)
            {
                final += item;
            }
            return final;
        }
    }
}
